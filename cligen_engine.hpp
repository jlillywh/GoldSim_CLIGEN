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

// Convert OLE automation date (days since Dec 30, 1899, as used by GoldSim DateTime)
// to Gregorian calendar (year, month 1..12, day 1..31, doy 1..366).
inline void ole_to_date(double ole_days, int& year, int& month, int& day, int& doy) {
    long jd = 2415019 + static_cast<long>(std::floor(ole_days + 1e-6));
    long l = jd + 68569;
    long n = (4 * l) / 146097;
    l = l - (146097 * n + 3) / 4;
    long i = (4000 * (l + 1)) / 1461001;
    l = l - (1461 * i) / 4 + 31;
    long j = (80 * l) / 2447;
    day = static_cast<int>(l - (2447 * j) / 80);
    l = j / 11;
    month = static_cast<int>(j + 2 - (12 * l));
    year = static_cast<int>(100 * (n - 49) + i + l);

    bool leap = is_leap_year(year);
    const int* days = leap ? DAYS_IN_MONTH_LEAP : DAYS_IN_MONTH;
    doy = day;
    for (int m = 0; m < month - 1; ++m) {
        doy += days[m];
    }
}

// ---------------------------------------------------------------------------
// Periodic (seasonal) cubic spline over the 12 monthly climate values.
//
// GoldSim now hands the generator a whole year of monthly values at once (one
// number per calendar month). This builds a smooth seasonal curve through those
// 12 anchor points and reads off the value for any given day, so a daily-step
// model no longer sees abrupt jumps at month boundaries. Each monthly value is
// treated as the climatological value at the CENTRE of its month, and the curve
// wraps around the calendar (December flows smoothly back into January), so the
// seasonal cycle has no seam. The interpolation is a true cubic spline with
// continuous first and second derivatives (periodic end conditions).
//
// Internally the 12 knots sit at integer positions 0..11 with period 12 and
// unit spacing; evaluate() takes a continuous "month-centre coordinate" u (see
// month_phase() below for the day-of-year mapping).
// ---------------------------------------------------------------------------
class MonthlyCubicSpline {
public:
    MonthlyCubicSpline() {
        for (int i = 0; i < 12; ++i) { y_[i] = 0.0; y2_[i] = 0.0; }
        valid_ = false;
    }

    // Supply the 12 monthly anchor values (Jan..Dec). Cheap no-op if unchanged,
    // so a model that holds parameters constant pays the solve cost only once,
    // while a slowly drifting climate simply rebuilds the curve as values move.
    void set(const double y[12]) {
        if (valid_) {
            bool same = true;
            for (int i = 0; i < 12; ++i) { if (y_[i] != y[i]) { same = false; break; } }
            if (same) return;
        }
        for (int i = 0; i < 12; ++i) y_[i] = y[i];
        solve();
        valid_ = true;
    }

    // Evaluate the seasonal curve at month-centre coordinate u (periodic, any
    // real value; u = m lands exactly on month m's anchor value).
    double evaluate(double u) const {
        u -= 12.0 * std::floor(u / 12.0);          // wrap into [0, 12)
        int i = static_cast<int>(std::floor(u));
        if (i < 0) i = 0;
        if (i > 11) i = 11;
        const int ip1 = (i + 1) % 12;
        const double t = u - i;                    // position within segment, [0,1)
        const double a = 1.0 - t, b = t;           // unit spacing (h = 1)
        return a * y_[i] + b * y_[ip1]
             + ((a * a * a - a) * y2_[i] + (b * b * b - b) * y2_[ip1]) / 6.0;
    }

private:
    double y_[12];   // monthly anchor values
    double y2_[12];  // second derivatives at the knots
    bool   valid_;

    // Second derivatives for a periodic cubic spline on 12 equally spaced knots.
    // With unit spacing the knot equations are, for every i (indices mod 12):
    //   y2[i-1] + 4 y2[i] + y2[i+1] = 6 (y[i+1] - 2 y[i] + y[i-1])
    // which is a symmetric circulant tridiagonal system (diagonal 4, off-diagonal
    // 1, plus the two wrap-around corners). Solved with the standard cyclic /
    // Sherman-Morrison reduction to two ordinary tridiagonal solves.
    void solve() {
        double rhs[12];
        for (int i = 0; i < 12; ++i) {
            const int im1 = (i + 11) % 12, ip1 = (i + 1) % 12;
            rhs[i] = 6.0 * (y_[ip1] - 2.0 * y_[i] + y_[im1]);
        }
        solve_cyclic(4.0, 1.0, rhs, y2_);
    }

    // x solves the N=12 circulant tridiagonal system with constant diagonal
    // `diag`, constant off-diagonals `off`, and corner entries `off` (periodic).
    static void solve_cyclic(double diag, double off, const double d[12], double x[12]) {
        const int N = 12;
        const double gamma = -diag;                // avoids cancellation in bb[0]
        double bb[12];
        for (int i = 0; i < N; ++i) bb[i] = diag;
        bb[0]     = diag - gamma;
        bb[N - 1] = diag - off * off / gamma;

        double y[12], q[12], u[12] = {0};
        solve_tridiag(bb, off, d, y);              // T y = d
        u[0] = gamma; u[N - 1] = off;
        solve_tridiag(bb, off, u, q);              // T q = u

        const double vy = y[0] + (off / gamma) * y[N - 1];
        const double vq = q[0] + (off / gamma) * q[N - 1];
        const double fact = vy / (1.0 + vq);
        for (int i = 0; i < N; ++i) x[i] = y[i] - fact * q[i];
    }

    // Thomas algorithm for an N=12 tridiagonal system with diagonal bb[i] and
    // constant sub/super-diagonal `off`.
    static void solve_tridiag(const double bb[12], double off, const double d[12], double x[12]) {
        const int N = 12;
        double cp[12], dp[12];
        cp[0] = off / bb[0];
        dp[0] = d[0] / bb[0];
        for (int i = 1; i < N; ++i) {
            const double m = bb[i] - off * cp[i - 1];
            cp[i] = off / m;
            dp[i] = (d[i] - off * dp[i - 1]) / m;
        }
        x[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i) x[i] = dp[i] - cp[i] * x[i + 1];
    }
};

// Map a day-of-year (1..365 or 366) to the spline's month-centre coordinate, placing
// each monthly anchor at the middle of its month and letting the ends wrap.
inline double month_phase(int doy, bool leap = false) {
    const int* days = leap ? DAYS_IN_MONTH_LEAP : DAYS_IN_MONTH;
    const int max_doy = leap ? 366 : 365;
    if (doy < 1) doy = 1;
    if (doy > max_doy) doy = max_doy;
    int cum = 0, m = 0;
    for (; m < 12; ++m) {
        if (doy <= cum + days[m]) break;
        cum += days[m];
    }
    if (m > 11) m = 11;
    const int day_in_month = doy - cum;                  // 1..days[m]
    const double f = (day_in_month - 0.5) / days[m];     // midpoint fraction 0..1
    return m + (f - 0.5);                                 // m at mid-month; eval() wraps
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

} // namespace cligen
