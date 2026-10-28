// goldsim_external.cpp
// ---------------------------------------------------------------------------
// Pure C++ GoldSim External element DLL for USDA CLIGEN weather generation.
//
// This replaces legacy Python bridges with a high-performance native DLL.
// There is NO Python runtime: GoldSim links directly to native C++ code.
//
// Exchange Pattern (STEP-BY-STEP SCALARS -- the only mode):
//   GoldSim evaluates the element on EVERY simulation timestep (DT = 1 day).
//   On each call the DLL generates exactly ONE day of weather and returns
//   three scalars: precipitation, daily maximum temperature and daily minimum
//   temperature. There is NO lookup table and no full-run buffer -- the model
//   and the generator march forward together, one day at a time.
//
//   Because the generator is driven step by step, it MUST remember its own
//   internal state between calls (the wet/dry chain, the day-to-day
//   temperature persistence residuals and the random stream). That state lives
//   in a single persistent session (see StepSession). Whenever GoldSim begins
//   a new Monte Carlo realization the session detects the change and reseeds,
//   so every realization is a fresh, independent weather sequence.
//
//   GoldSim may also feed time-varying CLIGEN parameters on each timestep
//   (e.g. a warming trend applied to the monthly means), and the returned day
//   simply reflects whatever parameters arrived for that step.
//
// Method IDs (GoldSim -> DLL):
//   0  XF_INITIALIZE       set up the session
//   1  XF_CALCULATE        generate ONE day, return 3 scalars
//   2  XF_REP_VERSION      report DLL version
//   3  XF_REP_ARGUMENTS    report input count and output count
//   99 XF_CLEANUP          tear down
//
// Argument Contract (inargs - exactly 13 scalar values per timestep):
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
// Outputs Contract (outargs - exactly 3 scalar values updated each timestep):
//   [0] Precip_Rate     Daily stochastic precipitation rate (mm/day or in/day, prcp / dt)
//   [1] Tmax            Daily maximum air temperature (deg C or deg F)
//   [2] Tmin            Daily minimum air temperature (deg C or deg F)
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

static constexpr int    NUM_STEP_INPUTS  = 13; // 4 controls + 5 precip + 4 temp scalars
static constexpr int    NUM_STEP_OUTPUTS = 3;  // Precip, Tmax, Tmin
static constexpr double DLL_VERSION      = 3.0;

// ---------------------------------------------------------------------------
// Persistent step-by-step session.
//
// This is what lets the generator be driven one day at a time: every piece of
// state that must survive from one GoldSim timestep to the next lives here --
// the wet/dry chain memory, the autoregressive temperature residuals and the
// random number stream. reset_realization() wipes that memory and reseeds, and
// is called automatically the first time a new realization index is seen so
// each realization is an independent sequence.
// ---------------------------------------------------------------------------
#include <cstdlib>
#include <string>
#include <cctype>

enum LogLevel {
    LOG_OFF   = 0,
    LOG_ERROR = 1,
    LOG_INFO  = 2,
    LOG_DEBUG = 3
};

inline LogLevel detect_log_level() {
    // 1. Check environment variable CLIGEN_LOG_LEVEL
    const char* env = std::getenv("CLIGEN_LOG_LEVEL");
    if (env) {
        std::string s(env);
        for (auto& c : s) c = std::toupper(static_cast<unsigned char>(c));
        if (s == "DEBUG" || s == "3") return LOG_DEBUG;
        if (s == "INFO"  || s == "2") return LOG_INFO;
        if (s == "ERROR" || s == "1") return LOG_ERROR;
        if (s == "OFF" || s == "NONE" || s == "0") return LOG_OFF;
    }
    // 2. Check for presence of cligen_debug.txt flag in working directory
    FILE* fdbg = fopen("cligen_debug.txt", "r");
    if (fdbg) {
        fclose(fdbg);
        return LOG_DEBUG;
    }
    // Default to LOG_ERROR so simulations run at maximum speed with zero I/O overhead
    return LOG_ERROR;
}

namespace {
    struct StepSession {
        bool initialized = false;
        long current_realization = -1;
        long step_count = 0;
        bool prev_wet = false;      // wet/dry Markov state carried from yesterday
        double res_tmax = 0.0;      // AR(1) max-temperature residual carried forward
        double res_tmin = 0.0;      // AR(1) min-temperature residual carried forward
        std::mt19937_64 rng;
        std::normal_distribution<double> norm{0.0, 1.0};
        std::uniform_real_distribution<double> uni{0.0, 1.0};
        LogLevel log_level = LOG_ERROR;
        FILE* log_file = nullptr;

        void open_log() {
            log_level = detect_log_level();
            if (log_level == LOG_OFF) return;

            // Only create/open the file upfront if INFO or DEBUG is requested
            if (log_level >= LOG_INFO && !log_file) {
                log_file = fopen("cligen_debug.log", "w");
                if (log_file) {
                    fprintf(log_file, "====================================================================================\n");
                    fprintf(log_file, " CligenWeather.dll v%.1f Log (Level: %s)\n", DLL_VERSION,
                        log_level == LOG_DEBUG ? "DEBUG" : (log_level == LOG_INFO ? "INFO" : "ERROR"));
                    fprintf(log_file, "====================================================================================\n");
                    fflush(log_file);
                }
            }
        }

        void log_error(const char* msg) {
            if (log_level < LOG_ERROR) return;
            if (!log_file) log_file = fopen("cligen_debug.log", "a");
            if (log_file) {
                fprintf(log_file, "[ERROR] Step %ld: %s\n", step_count, msg);
                fflush(log_file);
            }
        }

        void close_log() {
            if (log_file) {
                if (log_level >= LOG_INFO) {
                    fprintf(log_file, "[CligenWeather] Simulation finished. Total steps: %ld\n", step_count);
                }
                fclose(log_file);
                log_file = nullptr;
            }
        }

        void reset_realization(long real_id) {
            current_realization = real_id;
            uint64_t seed = (static_cast<uint64_t>(real_id) * 6364136223846793005ULL + 1442695040888963407ULL) ^ 0x9E3779B97F4A7C15ULL;
            if (seed == 0) seed = 189;
            rng.seed(seed);
            prev_wet = false;
            res_tmax = 0.0;
            res_tmin = 0.0;
            if (log_file && log_level >= LOG_INFO) {
                fprintf(log_file, "\n--- Realization %ld (Seed: %llu) ---\n", real_id, (unsigned long long)seed);
                fflush(log_file);
            }
        }
    };

    StepSession g_step_session;
}

// Single standard GoldSim External entry point.
GS_EXPORT void GS_CALL CligenWeather(int methodID, int* status,
                                     double* inargs, double* outargs) {
    if (status) *status = XF_SUCCESS;

    switch (methodID) {
    case XF_INITIALIZE: {
        g_step_session.open_log();
        g_step_session.initialized = true;
        g_step_session.current_realization = -1;
        g_step_session.step_count = 0;
        if (g_step_session.log_file && g_step_session.log_level >= LOG_INFO) {
            fprintf(g_step_session.log_file, "[CligenWeather] Session initialized.\n");
            fflush(g_step_session.log_file);
        }
        break;
    }

    case XF_REP_VERSION: {
        if (outargs) outargs[0] = DLL_VERSION;
        break;
    }

    case XF_REP_ARGUMENTS: {
        if (outargs) {
            outargs[0] = static_cast<double>(NUM_STEP_INPUTS);
            outargs[1] = static_cast<double>(NUM_STEP_OUTPUTS);
        }
        break;
    }

    case XF_CALCULATE: {
        if (!inargs || !outargs) {
            g_step_session.log_error("Null inargs or outargs pointer passed by GoldSim");
            if (status) *status = XF_FAILURE;
            return;
        }

        const double etime     = inargs[0];
        const double dt        = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        (void)etime; (void)dt;
        const long   real_id   = (inargs[2] >= 1.0) ? static_cast<long>(inargs[2]) : 1;
        const int    units     = (inargs[3] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        // Carry state forward between days; reseed only when a NEW realization
        // begins, giving each realization its own independent sequence.
        if (!g_step_session.initialized || g_step_session.current_realization != real_id) {
            g_step_session.open_log();
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

        // 1. Wet/Dry Markov Chain Occurrence (uses yesterday's state)
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
        outargs[0] = (dt > 0.0) ? (prcp / dt) : prcp; // Precip_Rate (mm/day or in/day)
        outargs[1] = tmax;                            // Tmax
        outargs[2] = tmin;                            // Tmin

        g_step_session.step_count++;
        if (g_step_session.log_file && g_step_session.log_level >= LOG_DEBUG) {
            if (g_step_session.step_count <= 100 || (g_step_session.step_count % 365 == 0)) {
                fprintf(g_step_session.log_file,
                    "Step %5ld | ETIME=%6.1f d | Real=%ld | In: PWW=%.3f PWD=%.3f MP=%.3f SDP=%.3f Skew=%.2f MTmax=%.1f SDTmax=%.2f MTmin=%.1f SDTmin=%.2f | Out: Wet=%d Precip_Rate=%7.4f %s Tmax=%5.2f Tmin=%5.2f\n",
                    g_step_session.step_count, etime, real_id,
                    pww, pwd, prcp_mean, prcp_sd, prcp_skew,
                    mean_tmax, sd_tmax, mean_tmin, sd_tmin,
                    (is_wet ? 1 : 0), outargs[0], (units == 1 ? "in/d" : "mm/d"),
                    outargs[1], outargs[2]);
                fflush(g_step_session.log_file);
            }
        }
        break;
    }

    case XF_CLEANUP: {
        g_step_session.close_log();
        g_step_session.initialized = false;
        break;
    }

    default:
        if (status) *status = XF_FAILURE;
        break;
    }
}

// Alias entry point so both "CligenWeather" and "cligen_weather_step" work in GoldSim
GS_EXPORT void GS_CALL cligen_weather_step(int methodID, int* status,
                                           double* inargs, double* outargs) {
    CligenWeather(methodID, status, inargs, outargs);
}
