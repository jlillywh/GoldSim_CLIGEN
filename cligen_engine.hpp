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
// The daily draw itself is draw_day(). GoldSim and CligenEngine::step both call it.
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

// One day's climate, already reduced to scalars (a month, or a spline
// evaluation for this day-of-year). units: 0 = metric (mm, deg C), 1 = US customary.
struct DayClimate {
    double pww = 0.0;
    double pwd = 0.0;
    double prcp_mean = 0.0;
    double prcp_sd = 0.0;
    double prcp_skew = 0.0;
    double mean_tmax = 0.0;
    double sd_tmax = 0.0;
    double mean_tmin = 0.0;
    double sd_tmin = 0.0;
    int units = 0;
};

// State that has to survive from one day to the next.
struct DayMemory {
    bool prev_wet = false;
    double res_tmax = 0.0;
    double res_tmin = 0.0;
};

struct DayWeather {
    bool is_wet = false;
    double prcp = 0.0;   // daily depth, same units as prcp_mean (not a rate)
    double tmax = 0.0;
    double tmin = 0.0;
};

// Single daily CLIGEN draw: Markov wet/dry, normal high and low with a
// mean-preserving wet-day shift, and Pearson Type III wet-day depth.
// Sanitizes the climate in place. Advances mem and the caller's random stream.
inline DayWeather draw_day(DayClimate& in, DayMemory& mem,
                           std::mt19937_64& rng,
                           std::normal_distribution<double>& norm,
                           std::uniform_real_distribution<double>& uni,
                           double r_cross = 0.65,
                           double rho_ar1 = 0.60) {
    if (in.pww <= 0.0 && in.pwd <= 0.0) {
        in.pww = 0.45;
        in.pwd = 0.20;
        in.prcp_mean = (in.units == 1) ? 0.25 : 6.0;
        in.prcp_sd   = (in.units == 1) ? 0.30 : 7.0;
        in.prcp_skew = 2.0;
        in.mean_tmax = (in.units == 1) ? 68.0 : 20.0;
        in.sd_tmax   = (in.units == 1) ? 7.2  : 4.0;
        in.mean_tmin = (in.units == 1) ? 50.0 : 10.0;
        in.sd_tmin   = (in.units == 1) ? 6.3  : 3.5;
    }
    const double min_t_sd = (in.units == 1) ? (MIN_TEMP_SD * 1.8) : MIN_TEMP_SD;
    const double min_p_sd = (in.units == 1) ? (MIN_PRCP_SD / 25.4) : MIN_PRCP_SD;
    in.sd_tmax   = std::max(min_t_sd, in.sd_tmax);
    in.sd_tmin   = std::max(min_t_sd, in.sd_tmin);
    in.prcp_sd   = std::max(min_p_sd, in.prcp_sd);
    in.prcp_mean = std::max((in.units == 1 ? 0.001 : 0.01), in.prcp_mean);
    in.prcp_skew = std::min(4.0, std::max(0.0, in.prcp_skew));

    const double p_trans = mem.prev_wet ? in.pww : in.pwd;
    const bool is_wet = uni(rng) < p_trans;
    mem.prev_wet = is_wet;

    const double z1 = norm(rng);
    const double z_ind = norm(rng);
    const double z2 = r_cross * z1 +
                      std::sqrt(std::max(0.0, 1.0 - r_cross * r_cross)) * z_ind;
    const double decay = std::sqrt(std::max(0.0, 1.0 - rho_ar1 * rho_ar1));
    mem.res_tmax = rho_ar1 * mem.res_tmax + decay * z1;
    mem.res_tmin = rho_ar1 * mem.res_tmin + decay * z2;

    // Wet days run cooler than dry days. The two shifts cancel, so the
    // unconditional mean of the high stays at mean_tmax.
    const double pw = std::max(0.01, std::min(0.99, in.pwd / std::max(0.01, 1.0 - in.pww + in.pwd)));
    const double delta_t = (in.units == 1) ? 3.6 : 2.0;
    const double t_shift = is_wet ? (-(1.0 - pw) * delta_t) : (pw * delta_t);
    double tmax = in.mean_tmax + t_shift + in.sd_tmax * mem.res_tmax;
    double tmin = in.mean_tmin + in.sd_tmin * mem.res_tmin;
    if (tmin > tmax) std::swap(tmin, tmax);

    double prcp = 0.0;
    if (is_wet) {
        const double z = norm(rng);
        if (std::fabs(in.prcp_skew) < 1e-4) {
            prcp = in.prcp_mean + in.prcp_sd * z;
        } else {
            const double inner = 1.0 + (in.prcp_skew * z) / 6.0 - (in.prcp_skew * in.prcp_skew) / 36.0;
            prcp = in.prcp_mean + in.prcp_sd * (2.0 / in.prcp_skew) * (inner * inner * inner - 1.0);
        }
        const double min_prcp = (in.units == 1) ? 0.004 : 0.1;
        const double round_factor = (in.units == 1) ? 10000.0 : 1000.0;
        prcp = std::max(min_prcp, std::round(prcp * round_factor) / round_factor);
    }

    DayWeather w;
    w.is_wet = is_wet;
    w.prcp = prcp;
    w.tmax = tmax;
    w.tmin = tmin;
    return w;
}

// ---------------------------------------------------------------------------
// USDA-ARS CLIGEN v5 (cligen.f) Sub-Daily Storm Disaggregation & Diurnal Temp
// ---------------------------------------------------------------------------

// Monthly mean dimensionless 30-min rainfall ratio wi(mo) = E[R_30 / R] (cligen.f R5MONB)
constexpr double CLIGEN_DEFAULT_WI[12] = {
    0.283, 0.283, 0.283, 0.304, 0.329, 0.393,
    0.486, 0.486, 0.393, 0.329, 0.283, 0.283
};

// 12-interval cumulative empirical distribution of dimensionless time-to-peak TIMPKD(0..12)
constexpr double CLIGEN_DEFAULT_TIMPKD[13] = {
    0.00, 0.08, 0.17, 0.25, 0.35, 0.45, 0.55,
    0.65, 0.74, 0.82, 0.89, 0.95, 1.00
};

// cligen.f subroutine DSTG(ai, k7): normalized Gamma PDF rejection sampler
// f(x) = x^6.28 * exp(6.28 * (1 - x)) with xn1 = 6.28 (B. Yu, 1999)
inline double cligen_dstg(double ai, std::mt19937_64& rng, std::uniform_real_distribution<double>& uni) {
    constexpr double xn1 = 6.28;
    const double xm = std::max(0.02, std::min(0.99, ai));
    for (int iter = 0; iter < 500; ++iter) {
        const double x = uni(rng);
        const double y = uni(rng);
        const double u = x / xm;
        if (u <= 0.0) continue;
        const double f = std::pow(u, xn1) * std::exp(xn1 * (1.0 - u));
        if (y <= f) return x;
    }
    return xm;
}

// cligen.f subroutine ALPHB(rain, wi, k7): dimensionless 30-min depth fraction r1 = alpha_0.5
// Bounded below by ab = 0.02083 and above by Arnold & Williams (1989) ajp = 1 - exp(-125 / R_mm).
inline double cligen_alphb(double rain_mm, double wi_mo, std::mt19937_64& rng, std::uniform_real_distribution<double>& uni) {
    constexpr double ab = 0.02083;
    double ajp = 1.0 - ab;
    if (rain_mm >= 25.4) {
        ajp = 1.0 - std::exp(-125.0 / rain_mm);
    }
    const double wi = std::max(ab, std::min(ajp, wi_mo));
    const double r = cligen_dstg(wi, rng, uni);
    return std::max(ab, std::min(ajp, r));
}

// cligen.f subroutine TIMEPK(timpkd, k10): dimensionless time-to-peak tp in (0.01, 0.99]
inline double cligen_timepk(std::mt19937_64& rng, std::uniform_real_distribution<double>& uni) {
    const double r = uni(rng);
    int i = 1;
    while (i < 12 && CLIGEN_DEFAULT_TIMPKD[i] < r) {
        ++i;
    }
    const double denom = std::max(1e-6, CLIGEN_DEFAULT_TIMPKD[i] - CLIGEN_DEFAULT_TIMPKD[i - 1]);
    double x = (static_cast<double>(i) - (CLIGEN_DEFAULT_TIMPKD[i] - r) / denom) / 12.0;
    if (x <= 0.01) x = 0.01;
    if (x > 0.99)  x = 0.99;
    return x;
}

// Solve WEPP/CLIGEN double-exponential parameter a > 0 from (1 - exp(-a)) / a = 1 / ip
inline double solve_cligen_double_exp_a(double ip) {
    const double target_ip = std::max(1.01, std::min(60.0, ip));
    double a = std::max(0.02, target_ip - 1.0);
    for (int iter = 0; iter < 16; ++iter) {
        const double ea = std::exp(-a);
        const double f = 1.0 - ea - a / target_ip;
        const double df = ea - 1.0 / target_ip;
        if (std::fabs(df) < 1e-12) break;
        const double step = f / df;
        a = std::max(1e-4, a - step);
        if (std::fabs(step) < 1e-8) break;
    }
    return a;
}

// Analytical cumulative normalized rainfall depth F(tau) in [0, 1] for WEPP/CLIGEN double-exponential storm
// Conserves 100% of storm depth: F(0) = 0, F(tp) = tp, F(1) = 1.0.
inline double cligen_double_exp_cum_fraction(double tau, double tp, double a) {
    if (tau <= 0.0) return 0.0;
    if (tau >= 1.0) return 1.0;
    const double exp_neg_a = std::exp(-a);
    const double denom = 1.0 - exp_neg_a;
    if (denom <= 1e-12) return tau;
    if (tau <= tp) {
        return tp * (std::exp(a * (tau / tp - 1.0)) - exp_neg_a) / denom;
    } else {
        return tp + (1.0 - tp) * (1.0 - std::exp(-a * (tau - tp) / (1.0 - tp))) / denom;
    }
}

struct SubdailyStorm {
    double dur   = 0.0; // storm duration in hours
    double tp    = 0.5; // dimensionless time-to-peak in (0.01, 0.99]
    double ip    = 1.5; // peak-to-average intensity ratio
    double a     = 1.0; // double-exponential shape parameter
    double start = 0.0; // storm start hour of day in [0, 24)
};

// Compute rainfall depth (in units of prcp) from storm `st` falling within [t0, t1] (hours from storm day's 00:00).
// For today's storm on [h0, h1], pass (h0, h1).
// For yesterday's storm continuing across midnight into today's [h0, h1], pass (h0 + 24.0, h1 + 24.0).
inline double eval_storm_interval_depth(double prcp, const SubdailyStorm& st, double t0, double t1) {
    if (prcp <= 0.0 || st.dur <= 0.0 || t1 <= t0) return 0.0;
    const double s_start = st.start;
    const double s_end   = s_start + st.dur;
    const double ov0     = std::max(t0, s_start);
    const double ov1     = std::min(t1, s_end);
    if (ov1 <= ov0) return 0.0;
    const double tau0 = (ov0 - s_start) / st.dur;
    const double tau1 = (ov1 - s_start) / st.dur;
    const double frac = std::max(0.0,
        cligen_double_exp_cum_fraction(tau1, st.tp, st.a) -
        cligen_double_exp_cum_fraction(tau0, st.tp, st.a));
    return prcp * frac;
}

// Draw CLIGEN sub-daily storm parameters (cligen.f DAY_GEN lines 3106-3138) using an isolated RNG stream
inline SubdailyStorm draw_cligen_storm(double prcp, int units, int month_idx, double dt_hours,
                                       std::mt19937_64& sub_rng,
                                       std::uniform_real_distribution<double>& uni) {
    SubdailyStorm st;
    if (prcp <= 0.0) return st;

    const int m = std::max(0, std::min(11, month_idx));
    const double rain_mm = (units == 1) ? (prcp * 25.4) : prcp;
    const double wi_mo = CLIGEN_DEFAULT_WI[m];

    // 1. Storm duration (hours): dur = 3.99 / (-2 * ln(1 - r1_dur))
    const double r1_dur = cligen_alphb(rain_mm, wi_mo, sub_rng, uni);
    const double rl_dur = std::max(0.04166, -2.0 * std::log(std::max(1e-6, 1.0 - r1_dur)));
    const double raw_dur = std::min(24.0, std::max(dt_hours, 3.99 / rl_dur));
    const int steps_per_day = std::max(1, static_cast<int>(std::round(24.0 / dt_hours)));
    const int storm_steps = std::max(1, std::min(steps_per_day, static_cast<int>(std::round(raw_dur / dt_hours))));
    st.dur = storm_steps * dt_hours;

    // 2. Dimensionless time-to-peak tp in (0.01, 0.99]
    st.tp = cligen_timepk(sub_rng, uni);

    // 3. Peak-to-average intensity ratio ip (cligen.f DAY_GEN lines 3128-3138)
    const double r1_peak = cligen_alphb(rain_mm, wi_mo, sub_rng, uni);
    const double rl_peak = -2.0 * std::log(std::max(1e-6, 1.0 - r1_peak));
    const double r5p = rain_mm * rl_peak;
    const double r_avg = rain_mm / std::max(0.0833, st.dur);
    const double tymax = 152.4 + 9.211 * rain_mm;
    const double xmav_max = std::max(1.01, tymax / std::max(0.01, r_avg));
    double ip = r5p / std::max(0.01, r_avg);
    if (ip > xmav_max) ip = xmav_max;
    if (ip < 1.01)     ip = 1.01;
    st.ip = ip;
    st.a  = solve_cligen_double_exp_a(ip);

    // 4. Storm start hour aligned to sub-daily step grid uniformly across all 24 hours [0, 24).
    //    Storms with start + dur > 24.0 continue across midnight (00:00) into the next day so
    //    all hours of the day (including midnight) have identical storm occurrence and peak probability.
    const int start_step = static_cast<int>(std::floor(uni(sub_rng) * steps_per_day));
    st.start = std::min(steps_per_day - 1, std::max(0, start_step)) * dt_hours;
    return st;
}

// Continuous 3-segment piecewise half-cosine diurnal temperature synthesis across midnight (00:00).
// Matches simulation_engine.js: Tmin occurs at 03:00 and Tmax occurs at 15:00 with C^1 continuity at 00:00.
inline double eval_diurnal_temp(double hour_of_day, double prev_tmax, double tmin, double tmax, double next_tmin) {
    constexpr double PI_VAL = 3.14159265358979323846;
    const double h = std::max(0.0, std::min(24.0, hour_of_day));
    if (h < 3.0) {
        // 00:00 to 03:00: cooling from yesterday's Tmax (at 15:00 yesterday) to today's Tmin (at 03:00)
        const double tau = (h + 9.0) / 12.0; // in [0.75, 1.0)
        return tmin + (prev_tmax - tmin) * 0.5 * (1.0 + std::cos(PI_VAL * tau));
    } else if (h < 15.0) {
        // 03:00 to 15:00: warming from today's Tmin (at 03:00) to today's Tmax (at 15:00)
        const double tau = (h - 3.0) / 12.0; // in [0.0, 1.0)
        return tmin + (tmax - tmin) * 0.5 * (1.0 - std::cos(PI_VAL * tau));
    } else {
        // 15:00 to 24:00: cooling from today's Tmax (at 15:00) toward tomorrow's Tmin (at 03:00 tomorrow)
        const double tau = (h - 15.0) / 12.0; // in [0.0, 0.75]
        return next_tmin + (tmax - next_tmin) * 0.5 * (1.0 + std::cos(PI_VAL * tau));
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

        DayClimate climate;
        climate.pww = p_.pww[m];
        climate.pwd = p_.pwd[m];
        climate.prcp_mean = p_.prcp_mean[m];
        climate.prcp_sd = p_.prcp_sd[m];
        climate.prcp_skew = p_.prcp_skew[m];
        climate.mean_tmax = p_.mean_tmax[m];
        climate.sd_tmax = p_.sd_tmax[m];
        climate.mean_tmin = p_.mean_tmin[m];
        climate.sd_tmin = p_.sd_tmin[m];
        climate.units = 0;

        DayMemory mem;
        mem.prev_wet = prev_wet_;
        mem.res_tmax = res_tmax_;
        mem.res_tmin = res_tmin_;
        DayWeather day = draw_day(climate, mem, rng_, norm_, uni_, p_.r_cross, p_.rho_ar1);
        prev_wet_ = mem.prev_wet;
        res_tmax_ = mem.res_tmax;
        res_tmin_ = mem.res_tmin;

        double tmax = day.tmax + warming_delta;
        double tmin = day.tmin + warming_delta;
        if (tmin > tmax) std::swap(tmin, tmax);

        // Storm shape is only needed for sub-daily disaggregation. The daily
        // GoldSim outputs do not use it, and it is drawn after the daily weather
        // so it does not change the precipitation or temperature sequence.
        double s_dur = 0.0, s_tp = 0.0, s_ip = 0.0, s_start = 0.0;
        if (day.is_wet) {
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
        r.prcp = day.prcp;
        r.is_wet = day.is_wet ? 1.0 : 0.0;
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
