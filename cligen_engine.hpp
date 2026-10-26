// cligen_engine.hpp
// ---------------------------------------------------------------------------
// Pure C++ implementation of USDA CLIGEN daily stochastic weather generation.
// Self-contained C++ engine with NO Python dependency, designed for GoldSim
// External Elements (see goldsim_external.cpp).
//
// CLIGEN Mathematical Formulation:
//   1. Precipitation Occurrence: First-order two-state Markov chain
//      using 12 monthly transition probabilities P(W|W) and P(W|D).
//   2. Precipitation Depth: 3-parameter Pearson Type III (skewed normal)
//      distribution on wet days via the Wilson-Hilferty transformation.
//   3. Daily Temperatures: 12 monthly maximum/minimum means and daily standard
//      deviations, conditioned on wet/dry occurrence, with cross-correlation
//      and lag-1 autoregressive persistence AR(1).
//   4. Sub-Daily Disaggregation: CLIGEN/WEPP double-exponential storm
//      hyetograph (exact daily mass balance) and Parton & Logan (1981)
//      diurnal temperature profile.
// ---------------------------------------------------------------------------
#pragma once

#include <cmath>
#include <cstdint>
#include <random>
#include <algorithm>
#include <vector>

namespace cligen {

// Minimum bounds on statistical parameters
constexpr double MIN_TEMP_SD  = 0.2;    // floor on daily temperature sd (deg)
constexpr double MIN_PRCP_SD  = 0.05;   // floor on daily precip sd (mm or in)

// Days-in-month tables (non-leap / leap).
constexpr int DAYS_IN_MONTH[12]      = {31,28,31,30,31,30,31,31,30,31,30,31};
constexpr int DAYS_IN_MONTH_LEAP[12] = {31,29,31,30,31,30,31,31,30,31,30,31};

inline bool is_leap_year(int y) {
    return (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0));
}

inline int doy_to_month(int doy, bool leap) {
    const int* days = leap ? DAYS_IN_MONTH_LEAP : DAYS_IN_MONTH;
    int cum = 0;
    for (int m = 0; m < 12; ++m) {
        cum += days[m];
        if (doy <= cum) return m;
    }
    return 11;
}

inline int date_to_doy(int year, int month, int day) {
    bool leap = is_leap_year(year);
    const int* days = leap ? DAYS_IN_MONTH_LEAP : DAYS_IN_MONTH;
    int m_idx = std::max(1, std::min(12, month)) - 1;
    int doy = 0;
    for (int m = 0; m < m_idx; ++m) doy += days[m];
    doy += std::max(1, std::min(days[m_idx], day));
    return doy;
}

inline void advance_calendar(int& year, int& month, int& day, int& doy) {
    bool leap = is_leap_year(year);
    const int* days = leap ? DAYS_IN_MONTH_LEAP : DAYS_IN_MONTH;
    day++;
    doy++;
    if (day > days[month - 1]) {
        day = 1;
        month++;
        if (month > 12) {
            month = 1;
            year++;
            doy = 1;
        }
    }
}

// Full monthly USDA CLIGEN parameter set.
struct Parameters {
    // 12 months x 5 columns: [P(W|W), P(W|D), Mean_P, SD_P, Skew_P]
    double pww[12]       = {0.45, 0.42, 0.48, 0.46, 0.40, 0.35, 0.32, 0.30, 0.34, 0.38, 0.42, 0.45};
    double pwd[12]       = {0.20, 0.18, 0.22, 0.20, 0.18, 0.14, 0.12, 0.11, 0.13, 0.16, 0.18, 0.20};
    double prcp_mean[12] = {6.5,  6.0,  7.5,  8.0,  7.0,  5.5,  4.5,  4.8,  5.2,  6.0,  6.8,  6.5};
    double prcp_sd[12]   = {7.2,  6.8,  8.5,  9.0,  8.2,  6.5,  5.5,  5.8,  6.2,  7.0,  7.8,  7.2};
    double prcp_skew[12] = {2.1,  2.0,  1.9,  1.8,  2.0,  2.2,  2.5,  2.4,  2.2,  2.0,  2.1,  2.1};

    // 12 months x 4 columns: [Mean_Tmax, SD_Tmax, Mean_Tmin, SD_Tmin]
    double mean_tmax[12] = {3.2,  6.1, 11.5, 16.2, 22.0, 27.5, 32.8, 31.5, 26.0, 18.5, 10.2,  4.0};
    double sd_tmax[12]   = {4.5,  4.6,  5.0,  5.2,  4.8,  4.5,  3.8,  3.9,  4.4,  4.8,  4.6,  4.4};
    double mean_tmin[12] = {-5.5,-3.8, -0.2,  3.5,  8.2, 13.0, 18.2, 17.5, 12.0,  5.8,  0.2, -4.8};
    double sd_tmin[12]   = {4.0,  3.8,  3.9,  4.0,  3.8,  3.6,  3.2,  3.3,  3.7,  3.8,  3.9,  4.0};

    // Temperature cross-correlation and persistence
    double r_cross = 0.65;  // Tmax/Tmin cross-correlation
    double rho_ar1 = 0.60;  // lag-1 daily persistence

    // Wet-day max temperature depression (deg C or deg F equivalent)
    double tmax_wet_depression[12] = {1.5, 1.8, 2.2, 2.5, 2.8, 3.0, 3.2, 3.0, 2.5, 2.0, 1.6, 1.5};

    // Sub-daily storm parameters (duration, peak intensity ratio, time-to-peak)
    double storm_dur_mean[12] = {6.0, 6.0, 6.0, 5.5, 5.0, 4.0, 3.0, 3.0, 4.0, 5.0, 6.0, 6.0};
    double storm_ip_mean[12]  = {3.5, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 7.0, 6.0, 5.0, 4.0, 3.5};
    double storm_tp_mean      = 0.35;
    double storm_tp_sd        = 0.18;
};

struct DayResult {
    double tmax;
    double tmin;
    double tmean;
    double dtr;
    double prcp;
    double is_wet;
    double is_frost;
    double is_hot;
    double mdd;

    // CLIGEN storm descriptors for sub-daily disaggregation (0 on dry days).
    double storm_dur   = 0.0;  // storm duration, hours
    double storm_tp    = 0.0;  // normalized time-to-peak (0..1)
    double storm_ip    = 0.0;  // peak-intensity ratio (peak/average)
    double storm_start = 0.0;  // storm start hour-of-day (0..24-dur)
};

class CligenEngine {
public:
    explicit CligenEngine(const Parameters& p = Parameters(), uint64_t seed = 1)
        : p_(p), rng_(seed), norm_(0.0, 1.0), uni_(0.0, 1.0) {
    }

    void reseed(uint64_t seed) {
        rng_.seed(seed);
        res_tmax_ = 0.0;
        res_tmin_ = 0.0;
        prev_wet_ = false;
    }

    // Advance one daily simulation step for calendar month (0..11).
    DayResult step(int month, double warming_delta = 0.0) {
        const int m = std::max(0, std::min(11, month));

        // 1. Wet/dry Markov chain occurrence
        double p_trans = prev_wet_ ? p_.pww[m] : p_.pwd[m];
        bool is_wet = uni_(rng_) < p_trans;
        prev_wet_ = is_wet;

        // 2. Daily Temperature Synthesis
        const double mu_max = is_wet ? (p_.mean_tmax[m] - p_.tmax_wet_depression[m]) : p_.mean_tmax[m];
        const double sd_max = std::max(MIN_TEMP_SD, p_.sd_tmax[m]);
        const double mu_min = p_.mean_tmin[m];
        const double sd_min = std::max(MIN_TEMP_SD, p_.sd_tmin[m]);

        // Cross-correlated bivariate normal variates
        const double z1 = norm_(rng_);
        const double z_ind = norm_(rng_);
        const double z2 = p_.r_cross * z1 +
                          std::sqrt(std::max(0.0, 1.0 - p_.r_cross * p_.r_cross)) * z_ind;

        // Lag-1 autoregressive persistence AR(1)
        const double decay = std::sqrt(std::max(0.0, 1.0 - p_.rho_ar1 * p_.rho_ar1));
        res_tmax_ = p_.rho_ar1 * res_tmax_ + decay * z1;
        res_tmin_ = p_.rho_ar1 * res_tmin_ + decay * z2;

        double tmax = mu_max + sd_max * res_tmax_ + warming_delta;
        double tmin = mu_min + sd_min * res_tmin_ + warming_delta;

        // Physical inversion swap guard
        if (tmin > tmax) std::swap(tmin, tmax);

        // 3. Precipitation Depth on Wet Days: Pearson Type III (Wilson-Hilferty)
        double prcp = 0.0;
        double s_dur = 0.0, s_tp = 0.0, s_ip = 0.0, s_start = 0.0;
        if (is_wet) {
            const double mean = std::max(0.01, p_.prcp_mean[m]);
            const double sd   = std::max(MIN_PRCP_SD, p_.prcp_sd[m]);
            const double g    = std::min(4.0, std::max(0.0, p_.prcp_skew[m]));
            const double z    = norm_(rng_);

            if (std::fabs(g) < 1e-4) {
                prcp = mean + sd * z;
            } else {
                const double inner = 1.0 + (g * z) / 6.0 - (g * g) / 36.0;
                prcp = mean + sd * (2.0 / g) * (inner * inner * inner - 1.0);
            }
            prcp = std::max(0.1, std::round(prcp * 1000.0) / 1000.0);

            // CLIGEN storm descriptors (duration, time-to-peak, peak intensity)
            const double dmean = std::max(0.25, p_.storm_dur_mean[m]);
            s_dur = -dmean * std::log(std::max(1e-9, uni_(rng_)));
            s_dur = std::min(24.0, std::max(0.25, s_dur));

            const double imean = std::max(0.1, p_.storm_ip_mean[m] - 1.0);
            s_ip = 1.0 - imean * std::log(std::max(1e-9, uni_(rng_)));
            s_ip = std::min(25.0, std::max(1.05, s_ip));

            s_tp = p_.storm_tp_mean + p_.storm_tp_sd * norm_(rng_);
            s_tp = std::min(0.95, std::max(0.05, s_tp));

            const double latest = std::max(0.0, 24.0 - s_dur);
            s_start = uni_(rng_) * latest;
        }

        DayResult r;
        r.tmax = tmax;
        r.tmin = tmin;
        r.tmean = 0.5 * (tmax + tmin);
        r.dtr = tmax - tmin;
        r.prcp = prcp;
        r.is_wet = is_wet ? 1.0 : 0.0;
        r.is_frost = (tmin <= 0.0) ? 1.0 : 0.0;
        r.is_hot = (tmax >= 32.22) ? 1.0 : 0.0;
        r.mdd = std::max(0.0, r.tmean);
        r.storm_dur = s_dur; r.storm_tp = s_tp; r.storm_ip = s_ip; r.storm_start = s_start;
        return r;
    }

private:
    Parameters p_;
    std::mt19937_64 rng_;
    std::normal_distribution<double> norm_;
    std::uniform_real_distribution<double> uni_;

    double res_tmax_ = 0.0;
    double res_tmin_ = 0.0;
    bool   prev_wet_ = false;
};

// Aliased for backwards compatibility in existing headers
using TemperatureEngine = CligenEngine;

// --- Full-run output as one lookup table PER WEATHER VARIABLE --------------
struct LookupTable {
    int table_id = 0;                 // which variable (see WxTableId)
    int table_dim = 2;
    std::vector<double> row_labels;   // day index, strictly increasing
    std::vector<double> col_labels;   // daily: variable codes; hourly: hour-of-day
    std::vector<double> data;         // row-major: data[row*ncols + col]
    int nrows = 0;
    int ncols = 0;
};

// Column codes used in the DAILY temperature/precip tables.
enum WxVar { WX_TMAX = 1, WX_TMIN = 2, WX_TMEAN = 3, WX_PRCP = 4 };
enum WxTableId { TBL_TEMPERATURE = 1, TBL_PRECIP = 2 };
enum StormMethod { STORM_CLIGEN_NATIVE = 0, STORM_METHOD_OF_FRAGMENTS = 1 };

// --- CLIGEN Sub-Daily Disaggregation --------------------------------------
inline double solve_decay_scale(double tp, double ip) {
    const double target = 1.0 / std::max(1.0001, ip);
    double lo = 1e-4, hi = 1e4;
    auto f = [&](double c) {
        return c * (2.0 - std::exp(-tp / c) - std::exp(-(1.0 - tp) / c));
    };
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

inline double storm_cum_depth(double t, double P, double dur,
                              double tp, double ip, double c) {
    if (t <= 0.0)   return 0.0;
    if (t >= dur)   return P;
    const double Tp = tp * dur;
    const double sc = c * dur;
    const double i_p = ip * P / dur;
    if (t <= Tp)
        return i_p * sc * (std::exp(-(Tp - t) / sc) - std::exp(-Tp / sc));
    const double Mr = i_p * sc * (1.0 - std::exp(-Tp / sc));
    return Mr + i_p * sc * (1.0 - std::exp(-(t - Tp) / sc));
}

inline double diurnal_temp(double hod, double tmin, double tmax) {
    constexpr double PI = 3.14159265358979323846;
    const double sunrise = 6.0, sunset = 20.0;
    const double daylen  = sunset - sunrise;
    const double lag     = 1.5;
    const double decay   = 2.2;
    auto day_t = [&](double h) {
        double arg = PI * (h - sunrise) / (daylen + 2.0 * lag);
        return tmin + (tmax - tmin) * std::sin(arg);
    };
    if (hod >= sunrise && hod <= sunset) return day_t(hod);
    const double tset     = day_t(sunset);
    const double nightlen = 24.0 - daylen;
    const double hs = (hod > sunset) ? (hod - sunset) : (hod + 24.0 - sunset);
    return tmin + (tset - tmin) * std::exp(-decay * hs / nightlen);
}

// Generate full simulation realization matching GoldSim start date and duration.
inline std::vector<LookupTable> generate_tables_cligen(
    const Parameters& p, uint64_t seed,
    int start_year, int start_month, int start_day,
    long duration_days, long steps_per_day, int storm_method,
    double warming_delta) {

    if (storm_method != STORM_CLIGEN_NATIVE) storm_method = STORM_CLIGEN_NATIVE;
    if (steps_per_day < 1) steps_per_day = 1;
    if (duration_days < 1) duration_days = 365;
    if (start_year < 1900) start_year = 2026;
    if (start_month < 1 || start_month > 12) start_month = 1;
    if (start_day < 1 || start_day > 31) start_day = 1;

    const bool hourly = (steps_per_day > 1);
    const double bin_hours = 24.0 / double(steps_per_day);

    CligenEngine eng(p, seed);

    LookupTable temp; temp.table_id = TBL_TEMPERATURE; temp.table_dim = 2;
    LookupTable rain; rain.table_id = TBL_PRECIP;      rain.table_dim = 2;

    if (!hourly) {
        temp.col_labels = { double(WX_TMAX), double(WX_TMIN) };
        rain.col_labels = { double(WX_PRCP) };
    } else {
        for (long h = 0; h < steps_per_day; ++h) {
            temp.col_labels.push_back(static_cast<double>(h));
            rain.col_labels.push_back(static_cast<double>(h));
        }
    }
    temp.ncols = static_cast<int>(temp.col_labels.size());
    rain.ncols = static_cast<int>(rain.col_labels.size());

    int cur_year  = start_year;
    int cur_month = start_month;
    int cur_day   = start_day;
    int cur_doy   = date_to_doy(cur_year, cur_month, cur_day);

    for (long d = 1; d <= duration_days; ++d) {
        DayResult r = eng.step(cur_month - 1, warming_delta);
        temp.row_labels.push_back(static_cast<double>(d));
        rain.row_labels.push_back(static_cast<double>(d));

        if (!hourly) {
            temp.data.push_back(r.tmax);
            temp.data.push_back(r.tmin);
            rain.data.push_back(r.prcp);
        } else {
            double c = 0.0;
            bool has_storm = (r.prcp > 0.0 && r.storm_dur > 0.0);
            if (has_storm) c = solve_decay_scale(r.storm_tp, r.storm_ip);
            for (long h = 0; h < steps_per_day; ++h) {
                double hod = bin_hours * (double(h) + 0.5);
                temp.data.push_back(diurnal_temp(hod, r.tmin, r.tmax));

                double depth = 0.0;
                if (has_storm) {
                    double a = std::max(bin_hours * double(h),       r.storm_start) - r.storm_start;
                    double b = std::min(bin_hours * double(h + 1),   r.storm_start + r.storm_dur) - r.storm_start;
                    if (b > a)
                        depth = storm_cum_depth(b, r.prcp, r.storm_dur, r.storm_tp, r.storm_ip, c)
                              - storm_cum_depth(a, r.prcp, r.storm_dur, r.storm_tp, r.storm_ip, c);
                }
                rain.data.push_back(std::max(0.0, depth));
            }
        }

        advance_calendar(cur_year, cur_month, cur_day, cur_doy);
    }

    temp.nrows = static_cast<int>(temp.row_labels.size());
    rain.nrows = static_cast<int>(rain.row_labels.size());
    return { temp, rain };
}

// Backward-compatible overload for integer simulation years.
inline std::vector<LookupTable> generate_tables(const Parameters& p, uint64_t seed,
                                int num_years, int start_year,
                                long steps_per_day, int storm_method,
                                double warming_delta) {
    long duration_days = std::max(1, num_years) * 365;
    return generate_tables_cligen(p, seed, start_year, 1, 1, duration_days,
                                  steps_per_day, storm_method, warming_delta);
}

// Flatten per-variable lookup tables into GoldSim contiguous outargs buffer.
inline long pack_tables(const std::vector<LookupTable>& tables,
                        double* outargs, long capacity) {
    long need = 1;
    for (const auto& t : tables)
        need += 4 + t.nrows + t.ncols + (long)t.nrows * t.ncols;
    if (!outargs || capacity < need) return -1;
    long k = 0;
    outargs[k++] = static_cast<double>(tables.size());
    for (const auto& t : tables) {
        outargs[k++] = t.table_id;
        outargs[k++] = t.table_dim;
        outargs[k++] = t.nrows;
        outargs[k++] = t.ncols;
        for (int i = 0; i < t.nrows; ++i) outargs[k++] = t.row_labels[i];
        for (int j = 0; j < t.ncols; ++j) outargs[k++] = t.col_labels[j];
        for (size_t i = 0; i < t.data.size(); ++i) outargs[k++] = t.data[i];
    }
    return k;
}

} // namespace cligen
