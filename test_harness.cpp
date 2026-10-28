// test_harness.cpp
// ---------------------------------------------------------------------------
// Standalone verification harness for pure C++ GoldSim Time Series generator.
// ---------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cassert>

// Forward declaration of the exported entry point
extern "C" void CligenWeather(int methodID, int* status,
                              double* inargs, double* outargs);

enum { XF_INITIALIZE = 0, XF_CALCULATE = 1, XF_REP_VERSION = 2, XF_REP_ARGUMENTS = 3, XF_CLEANUP = 99 };

int main(int argc, char** argv) {
    printf("===========================================================\n");
    printf(" Pure-C++ GoldSim External element (CLIGEN Time Series Mode)\n");
    printf("===========================================================\n");

    int status = 0;

    // 1. Check Reported Arguments
    double rep_out[8] = {0};
    CligenWeather(XF_REP_ARGUMENTS, &status, nullptr, rep_out);
    int num_inputs = static_cast<int>(rep_out[0]);
    int out_buffer_size = static_cast<int>(rep_out[1]);
    printf("[1] XF_REP_ARGUMENTS: inputs=%d (expected 113), buffer_size=%d (expected 500000)\n",
           num_inputs, out_buffer_size);
    assert(num_inputs == 113);
    assert(out_buffer_size >= 500000);

    // 2. Initialize
    CligenWeather(XF_INITIALIZE, &status, nullptr, nullptr);
    assert(status == 0);
    printf("[2] XF_INITIALIZE: status=%d (SUCCESS)\n", status);

    // 3. Prepare Inputs (113 values)
    // Salt Lake City test data (Oct 1, 2026 start date, 365 days duration, 1 day timestep)
    double inargs[113] = {0};
    inargs[0] = 46296.0; // StartTime (2026-10-01 OLE date)
    inargs[1] = 365.0;   // SimDuration (365 days)
    inargs[2] = 1.0;     // TimestepLength (1.0 day)
    inargs[3] = 1.0;     // Realization 1
    inargs[4] = 0.0;     // Units_Flag = 0 (Metric mm, deg C)

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

    // Pack the 9 monthly blocks into inargs[5..112]
    // 0: PWW, 1: PWD, 2: Mean_P, 3: SD_P, 4: Skew_P
    // 5: Mean_Tmax, 6: SD_Tmax, 7: Mean_Tmin, 8: SD_Tmin
    for (int m = 0; m < 12; ++m) {
        inargs[5   + m] = precip_table[m][0]; // PWW
        inargs[17  + m] = precip_table[m][1]; // PWD
        inargs[29  + m] = precip_table[m][2] * 25.4; // Mean_P (converted to mm)
        inargs[41  + m] = precip_table[m][3] * 25.4; // SD_P (converted to mm)
        inargs[53  + m] = precip_table[m][4]; // Skew_P
        inargs[65  + m] = temp_table[m][0];   // Mean_Tmax (deg C)
        inargs[77  + m] = temp_table[m][1];   // SD_Tmax (deg C)
        inargs[89  + m] = temp_table[m][2];   // Mean_Tmin (deg C)
        inargs[101 + m] = temp_table[m][3];   // SD_Tmin (deg C)
    }

    // Allocate output buffer
    std::vector<double> outargs(out_buffer_size, 0.0);

    // 4. Execute Calculation (Once for the entire realization!)
    CligenWeather(XF_CALCULATE, &status, inargs, outargs.data());
    printf("[3] XF_CALCULATE: status=%d (SUCCESS)\n", status);
    assert(status == 0);

    // 5. Verify the 3 Output Time Series Definitions
    const char* names[3] = {"Precip_Rate", "Tmax", "Tmin"};
    int expected_points = 366; // 365 days duration + 1 endpoint = 366 points
    int single_ts_size = 8 + 2 * expected_points;

    for (int s = 0; s < 3; ++s) {
        int offset = s * single_ts_size;
        double data_id   = outargs[offset + 0];
        double version   = outargs[offset + 1];
        double time_basis= outargs[offset + 2];
        double represents= outargs[offset + 3];
        double rows      = outargs[offset + 4];
        double cols      = outargs[offset + 5];
        double series    = outargs[offset + 6];
        double num_pts   = outargs[offset + 7];

        printf("  --- Output %d: %s ---\n", s + 1, names[s]);
        printf("    Format ID  = %.0f (expected 20)\n", data_id);
        printf("    Version    = %.0f (expected -3)\n", version);
        printf("    Time Basis = %.0f (expected 1: Calendar)\n", time_basis);
        printf("    Represents = %.0f (expected 0: Instantaneous)\n", represents);
        printf("    Rows/Cols  = %.0f / %.0f (scalar)\n", rows, cols);
        printf("    Histories  = %.0f (1)\n", series);
        printf("    Num Points = %.0f (expected %d)\n", num_pts, expected_points);

        assert(data_id == 20.0);
        assert(version == -3.0);
        assert(time_basis == 1.0);
        assert(represents == 0.0);
        assert(rows == 0.0 && cols == 0.0);
        assert(series == 1.0);
        assert(num_pts == static_cast<double>(expected_points));

        double t_first = outargs[offset + 8];
        double t_last  = outargs[offset + 8 + expected_points - 1];
        printf("    Time Span  = OLE %.1f to %.1f (delta = %.1f days)\n",
               t_first, t_last, t_last - t_first);
        assert(std::fabs(t_first - 46296.0) < 1e-4);
        assert(std::fabs(t_last - 46661.0) < 1e-4);

        // Value statistics
        double min_v = 1e9, max_v = -1e9, sum_v = 0.0;
        int val_offset = offset + 8 + expected_points;
        for (int i = 0; i < expected_points; ++i) {
            double v = outargs[val_offset + i];
            min_v = std::min(min_v, v);
            max_v = std::max(max_v, v);
            sum_v += v;
        }
        printf("    Values     : Min = %6.2f, Max = %6.2f, Avg = %6.2f, Total = %7.2f\n",
               min_v, max_v, sum_v / expected_points, sum_v);
    }

    // Check calendar temperature seasonality
    // Oct 1 (day 0) should be around 18-20 C
    // Jan (around day 100) should be coldest (around 0-5 C)
    // July (around day 290) should be hottest (around 30-35 C)
    int tmax_val_off = single_ts_size + 8 + expected_points;
    double tmax_oct = outargs[tmax_val_off + 0];   // Oct 1
    double tmax_jan = outargs[tmax_val_off + 105]; // Mid-Jan
    double tmax_jul = outargs[tmax_val_off + 290]; // Mid-July
    printf("\n[4] Seasonal Temperature Validation (Started Oct 1):\n");
    printf("    Oct 1 Tmax : %5.1f C\n", tmax_oct);
    printf("    Jan 15 Tmax: %5.1f C (winter trough)\n", tmax_jan);
    printf("    Jul 15 Tmax: %5.1f C (summer peak)\n", tmax_jul);

    assert(tmax_jul > tmax_oct);
    assert(tmax_oct > tmax_jan);

    // 6. Test Multi-Realization Seeding
    inargs[3] = 2.0; // Realization 2
    std::vector<double> outargs_real2(out_buffer_size, 0.0);
    CligenWeather(XF_CALCULATE, &status, inargs, outargs_real2.data());
    assert(status == 0);
    double r1_precip_sum = 0.0, r2_precip_sum = 0.0;
    int p_val_off = 8 + expected_points;
    for (int i = 0; i < expected_points; ++i) {
        r1_precip_sum += outargs[p_val_off + i];
        r2_precip_sum += outargs_real2[p_val_off + i];
    }
    printf("\n[5] Realization Variation:\n");
    printf("    Realization 1 Total Precip: %6.1f mm\n", r1_precip_sum);
    printf("    Realization 2 Total Precip: %6.1f mm\n", r2_precip_sum);
    assert(std::fabs(r1_precip_sum - r2_precip_sum) > 1.0); // Different stochastic realizations!

    // 7. Cleanup
    CligenWeather(XF_CLEANUP, &status, nullptr, nullptr);
    printf("\n[6] All tests PASSED successfully!\n");
    return 0;
}
