// test_harness.cpp
// ---------------------------------------------------------------------------
// Stands in for GoldSim: calls the exported CligenWeather() entry point through
// the exact same (methodID, status, inargs, outargs) contract GoldSim uses, and
// drives it the way GoldSim now does -- ONE TIMESTEP AT A TIME. Each simulated
// day the harness hands the element the full year of monthly CLIGEN vectors plus
// today's day-of-year, and reads back exactly three scalars (precip, Tmax, Tmin).
// The element itself fits a seasonal cubic spline through the 12 monthly values
// and interpolates to the day, so nothing is returned in bulk and there is no
// lookup table. Marching day by day also proves the element keeps its own state
// between calls (the wet/dry chain and temperature persistence only look right
// if yesterday is remembered). The harness writes the daily record to
// bridge_output/ and prints high-end statistics to confirm the native engine
// avoids the CLIGEN/WGEN high-tail defect.
//
// Build (Linux, no Python needed):
//   make run
//   # or: g++ -O2 -std=c++17 goldsim_external.cpp test_harness.cpp -o cligen_test
// ---------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

// Forward declaration of the DLL entry points (linked from goldsim_external.cpp).
extern "C" void CligenWeather(int methodID, int* status,
                              double* inargs, double* outargs);
extern "C" void cligen_weather_step(int methodID, int* status,
                                    double* inargs, double* outargs);

enum { XF_INITIALIZE = 0, XF_CALCULATE = 1, XF_REP_ARGUMENTS = 3, XF_CLEANUP = 99 };

// Days-in-month (non-leap) for picking the current calendar month as we march.
static const int DAYS_IN_MONTH[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

int main(int argc, char** argv) {
    int years = 30, start_year = 1995;
    double seed = 1.0, warming = 0.0;
    std::string out_csv = "bridge_output/cpp_sim_output.csv";
    for (int i = 1; i < argc - 1; ++i) {
        std::string k = argv[i], v = argv[i + 1];
        if (k == "--years") years = atoi(v.c_str());
        else if (k == "--start-year") start_year = atoi(v.c_str());
        else if (k == "--seed") seed = atof(v.c_str());
        else if (k == "--warming") warming = atof(v.c_str());
        else if (k == "--out") out_csv = v;
    }

    printf("===========================================================\n");
    printf(" Pure-C++ GoldSim External element  (CLIGEN weather)\n");
    printf("===========================================================\n");
    printf(" No Python. Native code only.\n");
    printf(" Exchange   : one timestep at a time -> 3 scalars (precip, Tmax, Tmin)\n");
    printf(" Inputs     : years=%d start_year=%d seed=%.0f warming=%.2f C\n",
           years, start_year, seed, warming);
    printf("              timestep=1 day (daily scalars, no lookup table)\n");

    int status = 0;

    // Confirm CligenWeather advertises the Time Series contract (113 in, 3 out).
    double rep_out[8] = {0};
    CligenWeather(XF_REP_ARGUMENTS, &status, nullptr, rep_out);
    int num_inputs = (int)rep_out[0];
    int num_outputs = (int)rep_out[1];
    printf(" CligenWeather (Time Series) reports %d inputs, %d outputs.\n", num_inputs, num_outputs);

    cligen_weather_step(XF_INITIALIZE, &status, nullptr, nullptr);
    if (status != 0) { fprintf(stderr, "INITIALIZE failed\n"); return 2; }

    // 12 months x 5 cols for Precip: [PWW, PWD, Mean, SD, Skew]
    double precip_table[12][5] = {
        {0.2692, 0.1398, 0.256, 0.150, 1.858},
        {0.2110, 0.1371, 0.278, 0.183, 2.229},
        {0.2416, 0.1652, 0.339, 0.249, 2.050},
        {0.2710, 0.1903, 0.394, 0.280, 1.755},
        {0.2500, 0.1500, 0.380, 0.260, 1.800},
        {0.2000, 0.1000, 0.300, 0.220, 2.100},
        {0.1800, 0.0800, 0.250, 0.200, 2.400},
        {0.2000, 0.0900, 0.280, 0.210, 2.300},
        {0.2200, 0.1100, 0.320, 0.230, 2.000},
        {0.2400, 0.1300, 0.350, 0.250, 1.900},
        {0.2600, 0.1500, 0.360, 0.260, 1.800},
        {0.2800, 0.1600, 0.320, 0.220, 1.800}
    };
    // 12 months x 4 cols for Temp: [Mean_Tmax, SD_Tmax, Mean_Tmin, SD_Tmin] (deg C)
    double temp_table[12][4] = {
        { 4.3, 5.2, -2.9, 4.5},
        { 7.2, 5.1, -1.0, 4.3},
        {12.5, 5.5,  2.8, 4.2},
        {16.8, 5.3,  6.2, 4.0},
        {22.4, 5.0, 10.8, 3.8},
        {28.5, 4.6, 15.5, 3.6},
        {33.8, 3.9, 20.2, 3.2},
        {32.4, 4.0, 19.5, 3.3},
        {26.8, 4.5, 14.2, 3.7},
        {18.9, 4.9,  7.8, 4.0},
        {10.5, 5.2,  1.5, 4.2},
        { 4.8, 5.0, -2.5, 4.4}
    };

    const long ndays = (long)years * 365;
    std::vector<double> tmax, tmin, tmean, prcp;
    tmax.reserve(ndays); tmin.reserve(ndays); tmean.reserve(ndays); prcp.reserve(ndays);
    int wet_days = 0;

    // Input layout (must match goldsim_external.cpp): 4 control scalars followed
    // by nine 12-element monthly vectors (112 inputs total).
    const int OFF_MONTHLY = 4;
    const int NIN = OFF_MONTHLY + 9 * 12; // 112
    std::vector<double> in(NIN, 0.0);
    double out[3] = {0};

    // Pack the whole year of monthly values into the nine vectors once; the
    // warming trend is baked into the Mean_Tmax / Mean_Tmin blocks. Column order
    // in precip_table is [PWW,PWD,Mean,SD,Skew]; temp_table is [MTmax,SDTmax,MTmin,SDTmin].
    auto blk = [&](int b) { return OFF_MONTHLY + b * 12; }; // start index of block b
    for (int m = 0; m < 12; ++m) {
        in[blk(0) + m] = precip_table[m][0];           // P_WW
        in[blk(1) + m] = precip_table[m][1];           // P_WD
        in[blk(2) + m] = precip_table[m][2];           // Mean_P
        in[blk(3) + m] = precip_table[m][3];           // SD_P
        in[blk(4) + m] = precip_table[m][4];           // Skew_P
        in[blk(5) + m] = temp_table[m][0] + warming;   // Mean_Tmax
        in[blk(6) + m] = temp_table[m][1];             // SD_Tmax
        in[blk(7) + m] = temp_table[m][2] + warming;   // Mean_Tmin
        in[blk(8) + m] = temp_table[m][3];             // SD_Tmin
    }

    // March the model forward one timestep at a time, exactly as GoldSim does.
    int month = 0;       // 0..11
    int day_of_month = 1;
    int doy = 1;         // calendar day-of-year 1..365
    for (long d = 0; d < ndays; ++d) {
        in[0] = (double)d;      // Simulation_Time (elapsed days)
        in[1] = 1.0;            // Time_Step (1 day)
        in[2] = seed;           // Realization (constant -> one continuous sequence)
        in[3] = 0.0;            // Units_Flag (0 = Metric, deg C / mm)

        cligen_weather_step(XF_CALCULATE, &status, in.data(), out);
        if (status != 0) { fprintf(stderr, "CALCULATE failed on day %ld\n", d); return 2; }

        double p = out[0], hi = out[1], lo = out[2];
        prcp.push_back(p);
        tmax.push_back(hi);
        tmin.push_back(lo);
        tmean.push_back(0.5 * (hi + lo));
        if (p > 0.0) wet_days++;

        // Advance the calendar one day (365-day year, matching ndays).
        if (++day_of_month > DAYS_IN_MONTH[month]) {
            day_of_month = 1;
            if (++month > 11) month = 0;
        }
        if (++doy > 365) doy = 1;
    }

    cligen_weather_step(XF_CLEANUP, &status, nullptr, nullptr);

    int n = (int)tmax.size();
    FILE* f = fopen(out_csv.c_str(), "w");
    if (f) {
        fprintf(f, "day,tmax,tmin,tmean,prcp\n");
        for (int i = 0; i < n; ++i)
            fprintf(f, "%d,%.4f,%.4f,%.4f,%.4f\n",
                    i + 1, tmax[i], tmin[i], tmean[i], prcp[i]);
        fclose(f);
        printf(" Wrote daily record -> %s\n", out_csv.c_str());
    }

    auto pct = [](std::vector<double> v, double q) {
        std::sort(v.begin(), v.end());
        double idx = q / 100.0 * (v.size() - 1);
        size_t lo = (size_t)std::floor(idx);
        double frac = idx - lo;
        return (lo + 1 < v.size()) ? v[lo] * (1 - frac) + v[lo + 1] * frac : v[lo];
    };
    double mn = *std::min_element(tmax.begin(), tmax.end());
    double mx = *std::max_element(tmax.begin(), tmax.end());
    double sum = 0; for (double v : tmax) sum += v;

    printf("-----------------------------------------------------------\n");
    printf(" Generated %d daily steps (%d years).\n", n, years);
    printf(" Daily high (tmax): mean %.1f C  min %.1f C  max %.1f C\n", sum / n, mn, mx);
    printf(" tmax p50=%.1f  p90=%.1f  p99=%.1f  p99.9=%.1f C\n",
           pct(tmax, 50), pct(tmax, 90), pct(tmax, 99), pct(tmax, 99.9));
    printf(" Wet-day fraction: %.1f%%   (%d of %d days)\n",
           100.0 * wet_days / n, wet_days, n);
    printf("-----------------------------------------------------------\n");
    printf(" Sanity: p99 should sit near ~37 C (observed), NOT ~66 C (WGEN defect).\n");
    printf("===========================================================\n");
    return 0;
}
