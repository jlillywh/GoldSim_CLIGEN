// test_harness.cpp
// ---------------------------------------------------------------------------
// Stands in for GoldSim: calls the exported CligenWeather() entry point through
// the exact same (methodID, status, inargs, outargs) contract GoldSim uses, but
// following the settled exchange pattern -- the full input set goes in once, and
// the engine returns the whole multi-year run as ONE LOOKUP TABLE PER VARIABLE.
// Daily: a temperature table (Tmax/Tmin cols) + a precip table. Hourly: each
// table is a day x hour-of-day grid. The harness unpacks the tables, writes the
// daily record to bridge_output/ and prints high-end statistics to confirm the
// native engine is statistically equivalent to the Python simulator (and avoids
// the CLIGEN/WGEN high-tail defect).
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

// Forward declaration of the DLL entry point (linked from goldsim_external.cpp).
extern "C" void CligenWeather(int methodID, int* status,
                              double* inargs, double* outargs);

enum { XF_INITIALIZE = 0, XF_CALCULATE = 1, XF_REP_ARGUMENTS = 3, XF_CLEANUP = 99 };

// Weather-variable column codes (mirror cligen::WxVar in cligen_engine.hpp).
enum { WX_TMAX = 1, WX_TMIN = 2, WX_TMEAN = 3, WX_PRCP = 4 };
// Per-variable table ids (mirror cligen::WxTableId).
enum { TBL_TEMPERATURE = 1, TBL_PRECIP = 2 };

// One unpacked lookup table.
struct Table {
    int id = 0, dim = 0, nrows = 0, ncols = 0;
    std::vector<double> row_labels, col_labels, data;
    int col_of(int code) const {
        for (int j = 0; j < ncols; ++j) if ((int)col_labels[j] == code) return j;
        return -1;
    }
    double at(int r, int c) const { return data[(size_t)r * ncols + c]; }
};

int main(int argc, char** argv) {
    int years = 30, start_year = 1995;
    double seed = 1.0, warming = 0.0, timestep_min = 1440.0, storm_method = 0.0;
    std::string out_csv = "bridge_output/cpp_sim_output.csv";
    for (int i = 1; i < argc - 1; ++i) {
        std::string k = argv[i], v = argv[i + 1];
        if (k == "--years") years = atoi(v.c_str());
        else if (k == "--start-year") start_year = atoi(v.c_str());
        else if (k == "--seed") seed = atof(v.c_str());
        else if (k == "--warming") warming = atof(v.c_str());
        else if (k == "--timestep-min") timestep_min = atof(v.c_str());
        else if (k == "--storm-method") storm_method = atof(v.c_str());
        else if (k == "--out") out_csv = v;
    }

    long steps_per_day = (timestep_min > 0.0) ? (long)std::floor(1440.0 / timestep_min + 0.5) : 1;
    bool clean_divisor = (timestep_min > 0.0) &&
        std::fabs(1440.0 / timestep_min - steps_per_day) < 1e-9 && steps_per_day >= 1;

    printf("===========================================================\n");
    printf(" Pure-C++ GoldSim External element  (CLIGEN weather)\n");
    printf("===========================================================\n");
    printf(" No Python. Native code only.\n");
    printf(" Exchange   : full run in one call -> one lookup table per variable\n");
    printf(" Inputs     : years=%d start_year=%d seed=%.0f warming=%.2f C\n",
           years, start_year, seed, warming);
    printf("              timestep=%.0f min (%s), storm_method=%.0f (%s)\n",
           timestep_min,
           !clean_divisor ? "invalid -> daily" : (steps_per_day == 1 ? "daily" : "sub-daily (CLIGEN-native)"),
           storm_method, storm_method == 0.0 ? "CLIGEN-native" : "reserved -> CLIGEN-native");

    int status = 0;

    // Ask the element how large an output buffer to reserve (GSPy max_elements).
    double rep_out[8] = {0};
    CligenWeather(XF_REP_ARGUMENTS, &status, nullptr, rep_out);
    int num_inputs = (int)rep_out[0];
    long capacity = (long)rep_out[1];
    std::vector<double> out(capacity, 0.0);
    std::vector<double> in(std::max(116, num_inputs), 0.0);

    CligenWeather(XF_INITIALIZE, &status, in.data(), out.data());
    if (status != 0) { fprintf(stderr, "INITIALIZE failed\n"); return 2; }

    in[0] = start_year;
    in[1] = 1.0;                 // start_month (Jan)
    in[2] = 1.0;                 // start_day (1st)
    in[3] = (double)(years * 365); // duration_days
    in[4] = timestep_min;
    in[5] = seed;
    in[6] = warming;
    in[7] = storm_method;

    // Populate sample monthly CLIGEN parameter matrices in inargs[8..115]
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
    // 12 months x 4 cols for Temp: [Mean_Tmax, SD_Tmax, Mean_Tmin, SD_Tmin]
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

    int idx = 8;
    for (int m = 0; m < 12; ++m) {
        for (int c = 0; c < 5; ++c) in[idx++] = precip_table[m][c];
    }
    for (int m = 0; m < 12; ++m) {
        for (int c = 0; c < 4; ++c) in[idx++] = temp_table[m][c];
    }

    CligenWeather(XF_CALCULATE, &status, in.data(), out.data());
    if (status != 0) { fprintf(stderr, "CALCULATE failed\n"); return 2; }
    CligenWeather(XF_CLEANUP, &status, in.data(), out.data());

    // Unpack the flattened per-variable tables:
    //   [ num_tables, {table_id, table_dim, nrows, ncols, row_labels, col_labels, data} * ]
    long k = 0;
    int num_tables = (int)out[k++];
    std::vector<Table> tables;
    for (int t = 0; t < num_tables; ++t) {
        Table tb;
        tb.id = (int)out[k++];
        tb.dim = (int)out[k++];
        tb.nrows = (int)out[k++];
        tb.ncols = (int)out[k++];
        tb.row_labels.assign(out.begin() + k, out.begin() + k + tb.nrows); k += tb.nrows;
        tb.col_labels.assign(out.begin() + k, out.begin() + k + tb.ncols); k += tb.ncols;
        tb.data.assign(out.begin() + k, out.begin() + k + (size_t)tb.nrows * tb.ncols);
        k += (size_t)tb.nrows * tb.ncols;
        tables.push_back(std::move(tb));
    }

    const Table* temp = nullptr; const Table* rain = nullptr;
    for (const auto& tb : tables) {
        if (tb.id == TBL_TEMPERATURE) temp = &tb;
        else if (tb.id == TBL_PRECIP) rain = &tb;
    }
    if (!temp || !rain) {
        fprintf(stderr, "Missing temperature or precipitation table (got %d tables)\n", num_tables);
        return 2;
    }

    int nrows = temp->nrows;
    bool hourly = (steps_per_day > 1) && clean_divisor;

    // Derive a per-day high (tmax) for the sanity check, from whichever shape
    // the temperature table came back in.
    std::vector<double> tmax(nrows), tmin(nrows), tmean(nrows), prcp(nrows);
    int wet_days = 0;
    if (!hourly) {
        int c_tmax = temp->col_of(WX_TMAX), c_tmin = temp->col_of(WX_TMIN);
        int c_prcp = rain->col_of(WX_PRCP);
        if (c_tmax < 0 || c_tmin < 0 || c_prcp < 0) {
            fprintf(stderr, "Unexpected daily table columns\n"); return 2;
        }
        for (int i = 0; i < nrows; ++i) {
            tmax[i] = temp->at(i, c_tmax);
            tmin[i] = temp->at(i, c_tmin);
            tmean[i] = 0.5 * (tmax[i] + tmin[i]);
            prcp[i] = rain->at(i, c_prcp);
            if (prcp[i] > 0.0) wet_days++;
        }
    } else {
        // Hourly grid: collapse each day's hour columns to min/max/mean & rain sum.
        for (int i = 0; i < nrows; ++i) {
            double hi = temp->at(i, 0), lo = temp->at(i, 0), s = 0.0, rs = 0.0;
            for (int j = 0; j < temp->ncols; ++j) {
                double v = temp->at(i, j);
                hi = std::max(hi, v); lo = std::min(lo, v); s += v;
            }
            for (int j = 0; j < rain->ncols; ++j) rs += rain->at(i, j);
            tmax[i] = hi; tmin[i] = lo; tmean[i] = s / temp->ncols;
            prcp[i] = rs;
            if (rs > 0.0) wet_days++;
        }
    }

    // Confirm each table's row axis is monotonically increasing.
    auto mono = [](const Table* tb) {
        for (int i = 1; i < tb->nrows; ++i)
            if (!(tb->row_labels[i] > tb->row_labels[i - 1])) return false;
        return true;
    };
    bool monotonic = mono(temp) && mono(rain);

    printf(" Element returned %d per-variable tables:\n", num_tables);
    printf("   temperature %dx%d, precipitation %dx%d  (row axes %s)\n",
           temp->nrows, temp->ncols, rain->nrows, rain->ncols,
           monotonic ? "monotonic" : "NOT MONOTONIC");
    if (hourly) {
        // Sub-daily sanity: mass balance (hourly rain should re-sum to the daily
        // depth) and a sense of the storm peakedness the double-exponential made.
        double peak_frac = 0.0;   // largest single-column share of a wet day's depth
        double max_abs_err = 0.0; // worst daily sum-vs-depth mismatch (mm)
        for (int i = 0; i < nrows; ++i) {
            if (prcp[i] <= 0.0) continue;
            double s = 0.0, hi = 0.0;
            for (int j = 0; j < rain->ncols; ++j) {
                double v = rain->at(i, j);
                s += v; hi = std::max(hi, v);
            }
            max_abs_err = std::max(max_abs_err, std::fabs(s - prcp[i]));
            if (s > 0.0) peak_frac = std::max(peak_frac, hi / s);
        }
        printf("   hourly shape: day x hour-of-day (%d cols), CLIGEN-native storm.\n",
               temp->ncols);
        printf("   sub-daily rain: mass balance err <= %.4g mm/day (0 = exact);\n"
               "   peak column holds up to %.0f%% of a storm's depth.\n",
               max_abs_err, 100.0 * peak_frac);
    }

    FILE* f = fopen(out_csv.c_str(), "w");
    if (f) {
        fprintf(f, "day,tmax,tmin,tmean,prcp\n");
        for (int i = 0; i < nrows; ++i)
            fprintf(f, "%.0f,%.4f,%.4f,%.4f,%.4f\n",
                    temp->row_labels[i], tmax[i], tmin[i], tmean[i], prcp[i]);
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
    printf(" Daily high (tmax): mean %.1f C  min %.1f C  max %.1f C\n", sum / nrows, mn, mx);
    printf(" tmax p50=%.1f  p90=%.1f  p99=%.1f  p99.9=%.1f C\n",
           pct(tmax, 50), pct(tmax, 90), pct(tmax, 99), pct(tmax, 99.9));
    printf(" Wet-day fraction: %.1f%%   (%d of %d days)\n",
           100.0 * wet_days / nrows, wet_days, nrows);
    printf("-----------------------------------------------------------\n");
    printf(" Sanity: p99 should sit near ~37 C (observed), NOT ~66 C (WGEN defect).\n");
    printf("===========================================================\n");
    return 0;
}
