// goldsim_external.cpp
// ---------------------------------------------------------------------------
// Pure C++ GoldSim External element DLL for USDA CLIGEN weather generation.
//
// This replaces legacy Python bridges with a high-performance native DLL.
// There is NO Python runtime: GoldSim links directly to native C++ code.
//
// Exchange Pattern:
//   GoldSim passes the simulation settings and the 12-month CLIGEN parameter
//   tables (defined in GoldSim Data elements) via inargs. At runtime, the DLL
//   synthesizes the requested realization for the exact simulation window
//   and returns self-describing 2D lookup tables to GoldSim.
//
// Method IDs (GoldSim -> DLL):
//   0  XF_INITIALIZE       set up the engine
//   1  XF_CALCULATE        generate the realization, return lookup tables
//   2  XF_REP_VERSION      report DLL version
//   3  XF_REP_ARGUMENTS    report input count and output capacity
//   99 XF_CLEANUP          tear down
//
// Argument Contract (inargs):
//   Control Parameters (8 values):
//     [0] start_year         calendar year (e.g. 2026)
//     [1] start_month        calendar start month (1..12, default 1)
//     [2] start_day          calendar start day (1..31, default 1)
//     [3] duration_days      simulation duration in days (e.g. 365, 3650)
//     [4] timestep_minutes   output resolution (1440 = daily, 60 = hourly)
//     [5] realization_seed   per-realization seed (GoldSim ~Realization index)
//     [6] warming_delta      climate-change temperature offset (deg C or deg F)
//     [7] storm_method       sub-daily rain shape (0 = CLIGEN double-exponential)
//
//   Monthly CLIGEN Parameters (108 values from GoldSim Data elements):
//     [8..67]   Precipitation Table (12 months x 5 columns):
//               For month m in 0..11:
//                 [8  + m*5 + 0] = P(W|W)   Wet after wet transition prob
//                 [8  + m*5 + 1] = P(W|D)   Wet after dry transition prob
//                 [8  + m*5 + 2] = Mean_P   Mean wet-day precip depth
//                 [8  + m*5 + 3] = SD_P     Std Dev of wet-day precip depth
//                 [8  + m*5 + 4] = Skew_P   Fisher-Pearson skewness coefficient
//
//     [68..115] Temperature Table (12 months x 4 columns):
//               For month m in 0..11:
//                 [68 + m*4 + 0] = Mean_Tmax  Monthly mean maximum daily temp
//                 [68 + m*4 + 1] = SD_Tmax    Daily Std Dev of maximum temp
//                 [68 + m*4 + 2] = Mean_Tmin  Monthly mean minimum daily temp
//                 [68 + m*4 + 3] = SD_Tmin    Daily Std Dev of minimum temp
//
// Outputs (outargs):
//   Self-describing 2D lookup tables (TBL_TEMPERATURE and TBL_PRECIP)
//   containing duration_days rows indexed 1..duration_days.
// ---------------------------------------------------------------------------
#include "cligen_engine.hpp"

#if defined(_WIN32)
  #define GS_EXPORT extern "C" __declspec(dllexport)
  #define GS_CALL   __stdcall
#else
  #define GS_EXPORT extern "C"
  #define GS_CALL
#endif

// GoldSim External method codes.
enum XFMethodID {
    XF_INITIALIZE       = 0,
    XF_CALCULATE        = 1,
    XF_REP_VERSION      = 2,
    XF_REP_ARGUMENTS    = 3,
    XF_CLEANUP          = 99
};

// GoldSim status codes (written back through *status).
enum XFStatus {
    XF_SUCCESS          = 0,
    XF_FAILURE          = 1,
    XF_FAILURE_WITH_MSG = -1
};

static constexpr int    NUM_CONTROL_INPUTS  = 8;
static constexpr int    NUM_PRECIP_INPUTS   = 60; // 12 x 5
static constexpr int    NUM_TEMP_INPUTS     = 48; // 12 x 4
static constexpr int    TOTAL_CLIGEN_INPUTS = NUM_CONTROL_INPUTS + NUM_PRECIP_INPUTS + NUM_TEMP_INPUTS; // 116
static constexpr double DLL_VERSION         = 2.1;

// Output buffer capacity GoldSim must reserve:
// ~100 years of daily rows (36,600 days), up to 24 sub-daily columns, for 2 tables.
static constexpr long   MAX_ROWS            = 36600;  // ~100 years of days
static constexpr long   MAX_COLS            = 24;     // hourly column ceiling
static constexpr long   MAX_OUTPUT_ELEMENTS =
    1 + 2L * (4 + MAX_ROWS + MAX_COLS + MAX_ROWS * MAX_COLS);

namespace {
    bool g_initialized = false;

    long steps_per_day_from_minutes(double timestep_minutes) {
        if (timestep_minutes <= 0.0) return 1;
        double steps = 1440.0 / timestep_minutes;
        double rounded = std::floor(steps + 0.5);
        if (std::fabs(steps - rounded) > 1e-9 || rounded < 1.0) return 1;
        return static_cast<long>(rounded);
    }
}

// Single standard GoldSim External entry point.
GS_EXPORT void GS_CALL CligenWeather(int methodID, int* status,
                                     double* inargs, double* outargs) {
    *status = XF_SUCCESS;

    switch (methodID) {
    case XF_INITIALIZE: {
        g_initialized = true;
        break;
    }

    case XF_REP_VERSION: {
        if (outargs) outargs[0] = DLL_VERSION;
        break;
    }

    case XF_REP_ARGUMENTS: {
        if (outargs) {
            outargs[0] = static_cast<double>(TOTAL_CLIGEN_INPUTS);
            outargs[1] = static_cast<double>(MAX_OUTPUT_ELEMENTS);
        }
        break;
    }

    case XF_CALCULATE: {
        if (!g_initialized) { *status = XF_FAILURE; return; }
        if (!inargs || !outargs) { *status = XF_FAILURE; return; }

        // Parse Control Parameters
        int start_year  = (inargs[0] >= 1900.0) ? static_cast<int>(inargs[0]) : 2026;
        int start_month = (inargs[1] >= 1.0 && inargs[1] <= 12.0) ? static_cast<int>(inargs[1]) : 1;
        int start_day   = (inargs[2] >= 1.0 && inargs[2] <= 31.0) ? static_cast<int>(inargs[2]) : 1;

        long duration_days = 0;
        if (inargs[3] > 0.0) {
            duration_days = static_cast<long>(inargs[3]);
        } else {
            // Backward compatibility fallback if years were passed in inargs[0]
            duration_days = std::max(1, static_cast<int>(inargs[0])) * 365;
        }

        const long     steps   = steps_per_day_from_minutes(inargs[4]);
        const uint64_t seed    = static_cast<uint64_t>(inargs[5] > 0.0 ? inargs[5] : 1.0);
        const double   warming = inargs[6];
        const int      storm   = static_cast<int>(inargs[7]);

        cligen::Parameters p;

        // Check if custom CLIGEN parameters were passed in from GoldSim Data elements
        bool has_custom_params = false;
        for (int i = NUM_CONTROL_INPUTS; i < TOTAL_CLIGEN_INPUTS; ++i) {
            if (inargs[i] != 0.0) {
                has_custom_params = true;
                break;
            }
        }

        if (has_custom_params) {
            int idx = NUM_CONTROL_INPUTS;
            for (int m = 0; m < 12; ++m) {
                p.pww[m]       = inargs[idx++];
                p.pwd[m]       = inargs[idx++];
                p.prcp_mean[m] = inargs[idx++];
                p.prcp_sd[m]   = inargs[idx++];
                p.prcp_skew[m] = inargs[idx++];
            }
            for (int m = 0; m < 12; ++m) {
                p.mean_tmax[m] = inargs[idx++];
                p.sd_tmax[m]   = inargs[idx++];
                p.mean_tmin[m] = inargs[idx++];
                p.sd_tmin[m]   = inargs[idx++];
            }
        }

        // Generate the realization matching GoldSim's start date and duration
        std::vector<cligen::LookupTable> tables = cligen::generate_tables_cligen(
            p, seed, start_year, start_month, start_day, duration_days, steps, storm, warming);

        long written = cligen::pack_tables(tables, outargs, MAX_OUTPUT_ELEMENTS);
        if (written < 0) {
            *status = XF_FAILURE;
            return;
        }
        break;
    }

    case XF_CLEANUP: {
        g_initialized = false;
        break;
    }

    default:
        *status = XF_FAILURE;
        break;
    }
}

// ---------------------------------------------------------------------------
// Step-by-Step Entry Point: cligen_weather_step
// ---------------------------------------------------------------------------
// Evaluated by GoldSim on EVERY simulation timestep (e.g. DT = 1 day).
// Enables GoldSim input parameters (e.g. PWW, PWD, Means) to be updated
// dynamically over simulation time (e.g. transient climate change trends).
//
// Arguments Contract (inargs - exactly 13 scalar values):
//   [0] Simulation_Time   Elapsed simulation time in days (GoldSim ETIME)
//   [1] Time_Step         Timestep size in days (GoldSim DT, default 1.0)
//   [2] Realization       Current Monte Carlo realization index (GoldSim Realization)
//   [3] Units_Flag        0 = Metric (mm, deg C), 1 = US Customary (in, deg F)
//   [4] P_WW              Wet-day prob given yesterday wet P(W|W) (-)
//   [5] P_WD              Wet-day prob given yesterday dry P(W|D) (-)
//   [6] Mean_P            Mean daily precipitation depth on wet days (mm or in)
//   [7] SD_P              Precipitation standard deviation (mm or in)
//   [8] Skew_P            Precipitation skewness coefficient (-)
//   [9] Mean_Tmax         Mean daily maximum air temperature (deg C or deg F)
//   [10] SD_Tmax          Maximum temperature standard deviation (Cdeg or Fdeg)
//   [11] Mean_Tmin        Mean daily minimum air temperature (deg C or deg F)
//   [12] SD_Tmin          Minimum temperature standard deviation (Cdeg or Fdeg)
//
// Outputs Contract (outargs - exactly 3 scalar values updated on each timestep):
//   [0] Precipitation   Daily stochastic rainfall depth (mm or in)
//   [1] Tmax            Daily maximum air temperature (deg C or deg F)
//   [2] Tmin            Daily minimum air temperature (deg C or deg F)
// ---------------------------------------------------------------------------
namespace {
    struct StepSession {
        bool initialized = false;
        long current_realization = -1;
        bool prev_wet = false;
        double res_tmax = 0.0;
        double res_tmin = 0.0;
        std::mt19937_64 rng;
        std::normal_distribution<double> norm{0.0, 1.0};
        std::uniform_real_distribution<double> uni{0.0, 1.0};

        void reset_realization(long real_id) {
            current_realization = real_id;
            uint64_t seed = (static_cast<uint64_t>(real_id) * 6364136223846793005ULL + 1442695040888963407ULL) ^ 0x9E3779B97F4A7C15ULL;
            if (seed == 0) seed = 189;
            rng.seed(seed);
            prev_wet = false;
            res_tmax = 0.0;
            res_tmin = 0.0;
        }
    };

    StepSession g_step_session;
}

GS_EXPORT void GS_CALL cligen_weather_step(int methodID, int* status,
                                           double* inargs, double* outargs) {
    if (status) *status = XF_SUCCESS;

    switch (methodID) {
    case XF_INITIALIZE: {
        g_step_session.initialized = true;
        g_step_session.current_realization = -1;
        break;
    }

    case XF_REP_VERSION: {
        if (outargs) outargs[0] = DLL_VERSION;
        break;
    }

    case XF_REP_ARGUMENTS: {
        if (outargs) {
            outargs[0] = 13.0; // 4 controls + 5 precip + 4 temp scalars
            outargs[1] = 3.0;  // exactly 3 outputs: Precip, Tmax, Tmin
        }
        break;
    }

    case XF_CALCULATE: {
        if (!inargs || !outargs) {
            if (status) *status = XF_FAILURE;
            return;
        }

        const double etime     = inargs[0];
        const double dt        = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        (void)etime; (void)dt;
        const long   real_id   = (inargs[2] >= 1.0) ? static_cast<long>(inargs[2]) : 1;
        const int    units     = (inargs[3] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        if (!g_step_session.initialized || g_step_session.current_realization != real_id) {
            g_step_session.initialized = true;
            g_step_session.reset_realization(real_id);
        }

        // 13 scalar inputs passed directly from GoldSim on this timestep
        double pww       = inargs[4];
        double pwd       = inargs[5];
        double prcp_mean = inargs[6];
        double prcp_sd   = inargs[7];
        double prcp_skew = inargs[8];

        double mean_tmax = inargs[9];
        double sd_tmax   = inargs[10];
        double mean_tmin = inargs[11];
        double sd_tmin   = inargs[12];

        // Sanity floors on variances
        if (pww <= 0.0 && pwd <= 0.0) {
            pww = 0.45; pwd = 0.20;
            prcp_mean = (units == 1) ? 0.25 : 6.0;
            prcp_sd   = (units == 1) ? 0.30 : 7.0;
            prcp_skew = 2.0;
            mean_tmax = (units == 1) ? 68.0 : 20.0;
            sd_tmax   = (units == 1) ? 7.2  : 4.0;
            mean_tmin = (units == 1) ? 50.0 : 10.0;
            sd_tmin   = (units == 1) ? 6.3  : 3.5;
        }
        const double min_t_sd = (units == 1) ? (cligen::MIN_TEMP_SD * 1.8) : cligen::MIN_TEMP_SD;
        const double min_p_sd = (units == 1) ? (cligen::MIN_PRCP_SD / 25.4) : cligen::MIN_PRCP_SD;
        sd_tmax = std::max(min_t_sd, sd_tmax);
        sd_tmin = std::max(min_t_sd, sd_tmin);
        prcp_sd = std::max(min_p_sd, prcp_sd);
        prcp_mean = std::max((units == 1 ? 0.001 : 0.01), prcp_mean);
        prcp_skew = std::min(4.0, std::max(0.0, prcp_skew));

        // 1. Wet/Dry Markov Chain Occurrence
        double p_trans = g_step_session.prev_wet ? pww : pwd;
        bool is_wet = g_step_session.uni(g_step_session.rng) < p_trans;
        g_step_session.prev_wet = is_wet;

        // 2. Temperature Synthesis (bivariate normal with cross-correlation and AR(1) persistence)
        constexpr double r_cross = 0.65;
        constexpr double rho_ar1 = 0.60;
        const double z1 = g_step_session.norm(g_step_session.rng);
        const double z_ind = g_step_session.norm(g_step_session.rng);
        const double z2 = r_cross * z1 + std::sqrt(1.0 - r_cross * r_cross) * z_ind;

        const double decay = std::sqrt(1.0 - rho_ar1 * rho_ar1);
        g_step_session.res_tmax = rho_ar1 * g_step_session.res_tmax + decay * z1;
        g_step_session.res_tmin = rho_ar1 * g_step_session.res_tmin + decay * z2;

        const double t_dep = is_wet ? (units == 1 ? 3.6 : 2.0) : 0.0;
        double tmax = mean_tmax - t_dep + sd_tmax * g_step_session.res_tmax;
        double tmin = mean_tmin + sd_tmin * g_step_session.res_tmin;
        if (tmin > tmax) std::swap(tmin, tmax);

        // 3. Precipitation Depth (Pearson Type III via Wilson-Hilferty transformation)
        double prcp = 0.0;
        if (is_wet) {
            const double z = g_step_session.norm(g_step_session.rng);
            if (std::fabs(prcp_skew) < 1e-4) {
                prcp = prcp_mean + prcp_sd * z;
            } else {
                const double inner = 1.0 + (prcp_skew * z) / 6.0 - (prcp_skew * prcp_skew) / 36.0;
                prcp = prcp_mean + prcp_sd * (2.0 / prcp_skew) * (inner * inner * inner - 1.0);
            }
            const double min_prcp = (units == 1) ? 0.004 : 0.1;
            const double round_factor = (units == 1) ? 10000.0 : 1000.0;
            prcp = std::max(min_prcp, std::round(prcp * round_factor) / round_factor);
        }

        // Return the 3 scalar outputs to GoldSim for the current timestep
        outargs[0] = prcp; // Precipitation
        outargs[1] = tmax; // Tmax
        outargs[2] = tmin; // Tmin
        break;
    }

    case XF_CLEANUP: {
        g_step_session.initialized = false;
        break;
    }

    default:
        if (status) *status = XF_FAILURE;
        break;
    }
}
