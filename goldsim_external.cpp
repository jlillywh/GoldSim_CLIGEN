// goldsim_external.cpp
// ---------------------------------------------------------------------------
// Pure C++ GoldSim External element DLL for USDA CLIGEN weather generation.
//
// This replaces legacy Python bridges with a high-performance native DLL.
// There is NO Python runtime: GoldSim links directly to native C++ code.
//
// Exchange Pattern (STEP-BY-STEP, MONTHLY VECTORS IN -- the only mode):
//   GoldSim evaluates the element on EVERY simulation timestep (DT = 1 day).
//   On each call the DLL generates exactly ONE day of weather and returns
//   three scalars: precipitation, daily maximum temperature and daily minimum
//   temperature. There is NO lookup table and no full-run buffer -- the model
//   and the generator march forward together, one day at a time.
//
//   The CLIGEN parameters are supplied as whole-year MONTHLY VECTORS: for each
//   of the nine climate parameters GoldSim passes 12 values (Jan..Dec). The DLL
//   fits a smooth periodic (seasonal) cubic spline through each set of 12
//   monthly anchors and evaluates it at the current day-of-year, so a daily
//   model sees a continuous seasonal signal rather than abrupt month-boundary
//   steps -- all interpolation happens inside the DLL, not in the GoldSim model.
//   GoldSim only needs to tell the element today's day-of-year (which it already
//   knows from its calendar); if that input is left at 0 the DLL derives the
//   day-of-year from elapsed simulation time assuming a 1-Jan start.
//
//   Because the generator is driven step by step, it MUST remember its own
//   internal state between calls (the wet/dry chain, the day-to-day
//   temperature persistence residuals and the random stream). That state lives
//   in a single persistent session (see StepSession). Whenever GoldSim begins
//   a new Monte Carlo realization the session detects the change and reseeds,
//   so every realization is a fresh, independent weather sequence.
//
//   A warming trend (or any climate drift) is applied by GoldSim evolving the
//   monthly vectors it passes; the DLL simply rebuilds the seasonal curve when
//   the incoming values change, so drift carries through automatically.
//
// Method IDs (GoldSim -> DLL):
//   0  XF_INITIALIZE       set up the session
//   1  XF_CALCULATE        generate ONE day, return 3 scalars
//   2  XF_REP_VERSION      report DLL version
//   3  XF_REP_ARGUMENTS    report input count and output count
//   99 XF_CLEANUP          tear down
//
// Argument Contract (inargs - exactly 112 values per timestep in CligenWeather / cligen_weather_arrays):
//   [0]  Simulation_Time   Elapsed simulation time in days (GoldSim ETIME)
//   [1]  Time_Step         Timestep size in days (GoldSim DT, default 1.0)
//   [2]  Realization       Current Monte Carlo realization index (GoldSim Realization)
//   [3]  Units_Flag        0 = Metric (mm, deg C), 1 = US Customary (in, deg F)
//   Then nine contiguous 12-element monthly vectors (Jan..Dec in each block):
//   [4..15]    P_WW[12]      Wet-day prob given yesterday wet P(W|W) (-)
//   [16..27]   P_WD[12]      Wet-day prob given yesterday dry P(W|D) (-)
//   [28..39]   Mean_P[12]    Mean daily precipitation depth on wet days (mm or in)
//   [40..51]   SD_P[12]      Precipitation standard deviation (mm or in)
//   [52..63]   Skew_P[12]    Precipitation skewness coefficient (-)
//   [64..75]   Mean_Tmax[12] Mean daily maximum air temperature (deg C or deg F)
//   [76..87]   SD_Tmax[12]   Maximum temperature standard deviation (Cdeg or Fdeg)
//   [88..99]   Mean_Tmin[12] Mean daily minimum air temperature (deg C or deg F)
//   [100..111] SD_Tmin[12]   Minimum temperature standard deviation (Cdeg or Fdeg)
//
// (Alternatively, cligen_weather_step accepts 13 scalar values for step-by-step models)
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

static constexpr int    NUM_STEP_OUTPUTS = 3;  // Precip, Tmax, Tmin
static constexpr int    NUM_TS_OUTPUTS   = 3;  // Precip_Rate, Tmax, Tmin Time Series
static constexpr double DLL_VERSION      = 5.0;

// GoldSim Time Series Definition protocol constants:
static constexpr int    EXTERNAL_TIMESERIES_DATA_ID = 20;
static constexpr int    EXTERNAL_TIMESERIES_VERSION = -3;
static int              g_ts_output_buffer_size     = 500000; // default buffer allocation in doubles

// Input layout for Time Series mode:
// 5 control scalars followed by nine 12-element monthly vectors (113 inputs total).
//   [0] Start_Time   (GoldSim StartTime / DateTime, OLE date double, e.g. 46296.0)
//   [1] Sim_Dur      (GoldSim SimDuration, days, e.g. 365.0)
//   [2] Time_Step    (GoldSim TimestepLength, days, e.g. 1.0)
//   [3] Realization  (GoldSim Realization, 1, 2, 3...)
//   [4] Units_Flag   (0 = Metric mm, deg C; 1 = US Customary in, deg F)
//   [5..112] 9 monthly vectors of 12 items
static constexpr int    MONTHLY_BLOCK    = 12;
static constexpr int    NUM_PARAM_BLOCKS = 9;   // PWW,PWD,MeanP,SDP,SkewP,MeanTmax,SDTmax,MeanTmin,SDTmin
static constexpr int    OFF_TS_MONTHLY   = 5;   // first monthly value sits at inargs[5]
static constexpr int    NUM_TS_INPUTS    = OFF_TS_MONTHLY + NUM_PARAM_BLOCKS * MONTHLY_BLOCK; // 113

// Input layout for legacy step-by-step array mode (112 inputs total):
static constexpr int    OFF_MONTHLY      = 4;   // first monthly value sits at inargs[4]
static constexpr int    NUM_ARRAY_INPUTS = OFF_MONTHLY + NUM_PARAM_BLOCKS * MONTHLY_BLOCK; // 112
static constexpr int    NUM_STEP_INPUTS  = 13;

// Start index of each monthly block within inargs.
enum ParamBlock {
    BLK_PWW = 0, BLK_PWD, BLK_MEANP, BLK_SDP, BLK_SKEWP,
    BLK_MTMAX, BLK_SDTMAX, BLK_MTMIN, BLK_SDTMIN
};
static inline int block_offset_step(int blk) { return OFF_MONTHLY + blk * MONTHLY_BLOCK; }
static inline int block_offset_ts(int blk)   { return OFF_TS_MONTHLY + blk * MONTHLY_BLOCK; }

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

#if defined(_WIN32)
#include <windows.h>
static std::string get_dll_directory() {
    char path[MAX_PATH] = {0};
    HMODULE hm = NULL;
    if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                           GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           (LPCSTR)&get_dll_directory, &hm)) {
        if (GetModuleFileNameA(hm, path, sizeof(path)) > 0) {
            char* slash = strrchr(path, '\\');
            if (slash) *slash = '\0';
            return std::string(path);
        }
    }
    return std::string(".");
}
#else
static std::string get_dll_directory() { return std::string("."); }
#endif

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
    // 2. Check for presence of cligen_debug.txt flag in working directory or DLL directory
    FILE* fdbg = fopen("cligen_debug.txt", "r");
    if (fdbg) {
        fclose(fdbg);
        return LOG_DEBUG;
    }
    std::string dll_dbg = get_dll_directory() + "\\cligen_debug.txt";
    fdbg = fopen(dll_dbg.c_str(), "r");
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
        // One seasonal cubic spline per CLIGEN parameter, rebuilt only when the
        // incoming monthly vectors change (so steady runs pay the solve once).
        cligen::MonthlyCubicSpline param_spline[NUM_PARAM_BLOCKS];
        std::mt19937_64 rng;
        std::normal_distribution<double> norm{0.0, 1.0};
        std::uniform_real_distribution<double> uni{0.0, 1.0};
        LogLevel log_level = LOG_ERROR;
        FILE* log_file = nullptr;
        std::string log_filepath;

        void open_log() {
            log_level = detect_log_level();
            if (log_level == LOG_OFF) return;

            if (log_filepath.empty()) {
                std::string dll_dir = get_dll_directory();
                log_filepath = dll_dir + "\\cligen_debug.log";
            }

            // Only create/open the file upfront if INFO or DEBUG is requested
            if (log_level >= LOG_INFO && !log_file) {
                log_file = fopen(log_filepath.c_str(), "w");
                if (!log_file) log_file = fopen("cligen_debug.log", "w");
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
            if (!log_file) {
                if (!log_filepath.empty()) log_file = fopen(log_filepath.c_str(), "a");
                if (!log_file) log_file = fopen("cligen_debug.log", "a");
            }
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

// Common execution function supporting both 112-input vector spline mode and 13-input scalar mode.
static void execute_cligen(bool array_mode, int methodID, int* status,
                           double* inargs, double* outargs) {
    if (status) *status = XF_SUCCESS;

    switch (methodID) {
    case XF_INITIALIZE: {
        g_step_session.open_log();
        g_step_session.initialized = true;
        g_step_session.current_realization = -1;
        g_step_session.step_count = 0;
        if (g_step_session.log_file && g_step_session.log_level >= LOG_INFO) {
            fprintf(g_step_session.log_file, "[CligenWeather] Session initialized (%s mode).\n",
                array_mode ? "112-vector spline" : "13-scalar step");
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
            outargs[0] = static_cast<double>(array_mode ? NUM_ARRAY_INPUTS : NUM_STEP_INPUTS);
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

        const double time_val  = inargs[0]; // GoldSim DateTime run property (preferred) or ETIME
        const double dt        = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        const long   real_id   = (inargs[2] >= 1.0) ? static_cast<long>(inargs[2]) : 1;
        const int    units     = (inargs[3] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        // Carry state forward between days; reseed only when a NEW realization begins
        if (!g_step_session.initialized || g_step_session.current_realization != real_id) {
            g_step_session.open_log();
            g_step_session.initialized = true;
            g_step_session.reset_realization(real_id);
        }

        double pww, pwd, prcp_mean, prcp_sd, prcp_skew;
        double mean_tmax, sd_tmax, mean_tmin, sd_tmin;

        if (array_mode) {
            int doy = 1;
            bool leap = false;
            if (time_val >= 10000.0) {
                // GoldSim DateTime run property (OLE automation date, days since 1899-12-30)
                int year = 0, month = 0, day = 0;
                cligen::ole_to_date(time_val, year, month, day, doy);
                leap = cligen::is_leap_year(year);
            } else {
                // Elapsed simulation time ETIME starting from day 0 (Jan 1)
                long day_idx = static_cast<long>(std::floor(time_val / dt + 0.5));
                if (day_idx < 0) day_idx = 0;
                doy = static_cast<int>((day_idx % 365) + 1);
            }

            // Refresh the seasonal splines from the monthly vectors and evaluate today's value
            const double u = cligen::month_phase(doy, leap);
            double param[NUM_PARAM_BLOCKS];
            for (int blk = 0; blk < NUM_PARAM_BLOCKS; ++blk) {
                g_step_session.param_spline[blk].set(&inargs[block_offset_step(blk)]);
                param[blk] = g_step_session.param_spline[blk].evaluate(u);
            }

            pww       = std::max(0.01, std::min(0.99, param[BLK_PWW]));
            pwd       = std::max(0.01, std::min(0.99, param[BLK_PWD]));
            prcp_mean = param[BLK_MEANP];
            prcp_sd   = param[BLK_SDP];
            prcp_skew = std::max(0.0, std::min(4.0, param[BLK_SKEWP]));

            mean_tmax = param[BLK_MTMAX];
            sd_tmax   = param[BLK_SDTMAX];
            mean_tmin = param[BLK_MTMIN];
            sd_tmin   = param[BLK_SDTMIN];
        } else {
            // Direct daily scalar inputs
            pww       = inargs[4];
            pwd       = inargs[5];
            prcp_mean = inargs[6];
            prcp_sd   = inargs[7];
            prcp_skew = inargs[8];

            mean_tmax = inargs[9];
            sd_tmax   = inargs[10];
            mean_tmin = inargs[11];
            sd_tmin   = inargs[12];
        }

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
        sd_tmax   = std::max(min_t_sd, sd_tmax);
        sd_tmin   = std::max(min_t_sd, sd_tmin);
        prcp_sd   = std::max(min_p_sd, prcp_sd);
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

        // Wet/dry temperature conditioning: wet days are cooled relative to dry days
        // while preserving the unconditional monthly mean (mean_tmax).
        // Steady-state wet probability from the two-state Markov chain:
        const double pw = std::max(0.01, std::min(0.99, pwd / std::max(0.01, 1.0 - pww + pwd)));
        const double delta_t = (units == 1) ? 3.6 : 2.0;
        // Dry days are raised by pw * delta_t; wet days are lowered by (1 - pw) * delta_t.
        // Expected value: pw * (mean - (1-pw)*delta) + (1-pw) * (mean + pw*delta) == mean_tmax.
        const double t_shift = is_wet ? (-(1.0 - pw) * delta_t) : (pw * delta_t);
        double tmax = mean_tmax + t_shift + sd_tmax * g_step_session.res_tmax;
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
                    "Step %5ld | Time=%6.1f | Real=%ld | In: PWW=%.3f PWD=%.3f MP=%.3f SDP=%.3f Skew=%.2f MTmax=%.1f SDTmax=%.2f MTmin=%.1f SDTmin=%.2f | Out: Wet=%d Precip_Rate=%7.4f %s Tmax=%5.2f Tmin=%5.2f\n",
                    g_step_session.step_count, time_val, real_id,
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

// ---------------------------------------------------------------------------
// Time Series Generation Mode (GoldSim Time Series Definition protocol):
// ---------------------------------------------------------------------------
// GoldSim evaluates this external function ONCE per realization.
// It generates the full synthetic weather sequence from Start_Time to
// (Start_Time + Sim_Dur) and returns three complete Time Series Definitions:
//   [0] Precip_Rate (mm/d or in/d, instantaneous values)
//   [1] Tmax        (deg C or deg F, instantaneous values)
//   [2] Tmin        (deg C or deg F, instantaneous values)
//
// Format for each scalar time series returned to GoldSim:
//   dest[0] = 20.0; // EXTERNAL_TIMESERIES_DATA_ID = 20
//   dest[1] = -3.0; // EXTERNAL_TIMESERIES_VERSION = -3
//   dest[2] = 1.0;  // 1 = Calendar time basis (OLE Julian date)
//   dest[3] = 0.0;  // 0 = Instantaneous values (interpolated by GoldSim)
//   dest[4] = 0.0;  // 0 rows (scalar)
//   dest[5] = 0.0;  // 0 cols (scalar)
//   dest[6] = 1.0;  // 1 series
//   dest[7] = num_points;
//   dest[8 .. 8+num_points-1] = timestamps (calendar dates: start_time + k*dt)
//   dest[8+num_points .. 8+2*num_points-1] = values
// ---------------------------------------------------------------------------
#include <cstring>
#include <numeric>

static inline void write_time_series_record(
    double* dest,
    int time_basis,     // 1 = Calendar, 0 = Elapsed
    int data_type,      // 0 = Instantaneous, 1 = Constant next
    int num_points,
    const double* timestamps,
    const double* values)
{
    dest[0] = static_cast<double>(EXTERNAL_TIMESERIES_DATA_ID); // 20.0
    dest[1] = static_cast<double>(EXTERNAL_TIMESERIES_VERSION); // -3.0
    dest[2] = static_cast<double>(time_basis); // 1.0
    dest[3] = static_cast<double>(data_type);  // 0.0
    dest[4] = 0.0;
    dest[5] = 0.0;
    dest[6] = 1.0;
    dest[7] = static_cast<double>(num_points);
    std::memcpy(&dest[8], timestamps, num_points * sizeof(double));
    std::memcpy(&dest[8 + num_points], values, num_points * sizeof(double));
}

static void execute_cligen_timeseries(int methodID, int* status,
                                      double* inargs, double* outargs) {
    if (status) *status = XF_SUCCESS;

    switch (methodID) {
    case XF_INITIALIZE: {
        g_step_session.open_log();
        g_step_session.initialized = true;
        g_step_session.current_realization = -1;
        g_step_session.step_count = 0;
        if (g_step_session.log_file && g_step_session.log_level >= LOG_INFO) {
            fprintf(g_step_session.log_file, "[CligenWeather] Time Series Generator initialized.\n");
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
            outargs[0] = static_cast<double>(NUM_TS_INPUTS); // 113 inputs
            outargs[1] = static_cast<double>(g_ts_output_buffer_size); // 500,000 doubles initial buffer
        }
        break;
    }

    case XF_CALCULATE: {
        if (!inargs || !outargs) {
            g_step_session.log_error("Null inargs or outargs pointer passed by GoldSim");
            if (status) *status = XF_FAILURE;
            return;
        }

        const double start_time = inargs[0]; // GoldSim StartTime / DateTime run property
        const double sim_dur    = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        const double dt         = (inargs[2] > 0.0) ? inargs[2] : 1.0;
        const long   real_id    = (inargs[3] >= 1.0) ? static_cast<long>(std::round(inargs[3])) : 1;
        const int    units      = (inargs[4] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        int num_steps = static_cast<int>(std::round(sim_dur / dt));
        if (num_steps < 1) num_steps = 1;
        int num_points = num_steps + 1; // from start_time to start_time + sim_dur inclusive

        int single_ts_size = 8 + 2 * num_points;
        int total_req_size = NUM_TS_OUTPUTS * single_ts_size;

        if (total_req_size > g_ts_output_buffer_size) {
            outargs[0] = static_cast<double>(total_req_size);
            g_ts_output_buffer_size = total_req_size;
            if (status) *status = -2; // XF_INCREASE_MEMORY
            return;
        }

        // Initialize realization RNG
        g_step_session.open_log();
        g_step_session.initialized = true;
        g_step_session.reset_realization(real_id);

        // Build seasonal splines from the 9 monthly vectors
        cligen::MonthlyCubicSpline param_spline[NUM_PARAM_BLOCKS];
        for (int blk = 0; blk < NUM_PARAM_BLOCKS; ++blk) {
            param_spline[blk].set(&inargs[block_offset_ts(blk)]);
        }

        std::vector<double> ts_times(num_points);
        std::vector<double> ts_prcp(num_points);
        std::vector<double> ts_tmax(num_points);
        std::vector<double> ts_tmin(num_points);

        bool prev_wet = false;
        double res_tmax = 0.0;
        double res_tmin = 0.0;

        for (int k = 0; k < num_points; ++k) {
            double cur_time = start_time + k * dt;
            ts_times[k] = cur_time;

            int doy = 1;
            bool leap = false;
            if (cur_time >= 10000.0) {
                int year = 0, month = 0, day = 0;
                cligen::ole_to_date(cur_time, year, month, day, doy);
                leap = cligen::is_leap_year(year);
            } else {
                long day_idx = static_cast<long>(std::floor(cur_time / dt + 0.5));
                if (day_idx < 0) day_idx = 0;
                doy = static_cast<int>((day_idx % 365) + 1);
            }

            const double u = cligen::month_phase(doy, leap);
            double param[NUM_PARAM_BLOCKS];
            for (int blk = 0; blk < NUM_PARAM_BLOCKS; ++blk) {
                param[blk] = param_spline[blk].evaluate(u);
            }

            double pww       = std::max(0.01, std::min(0.99, param[BLK_PWW]));
            double pwd       = std::max(0.01, std::min(0.99, param[BLK_PWD]));
            double prcp_mean = param[BLK_MEANP];
            double prcp_sd   = param[BLK_SDP];
            double prcp_skew = std::max(0.0, std::min(4.0, param[BLK_SKEWP]));

            double mean_tmax = param[BLK_MTMAX];
            double sd_tmax   = param[BLK_SDTMAX];
            double mean_tmin = param[BLK_MTMIN];
            double sd_tmin   = param[BLK_SDTMIN];

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
            sd_tmax   = std::max(min_t_sd, sd_tmax);
            sd_tmin   = std::max(min_t_sd, sd_tmin);
            prcp_sd   = std::max(min_p_sd, prcp_sd);
            prcp_mean = std::max((units == 1 ? 0.001 : 0.01), prcp_mean);
            prcp_skew = std::min(4.0, std::max(0.0, prcp_skew));

            // 1. Wet/Dry Markov chain occurrence
            double p_trans = prev_wet ? pww : pwd;
            bool is_wet = g_step_session.uni(g_step_session.rng) < p_trans;
            prev_wet = is_wet;

            // 2. Temperature residuals (Matalas AR(1))
            constexpr double r_cross = 0.65;
            constexpr double rho_ar1 = 0.60;
            const double z1 = g_step_session.norm(g_step_session.rng);
            const double z_ind = g_step_session.norm(g_step_session.rng);
            const double z2 = r_cross * z1 + std::sqrt(1.0 - r_cross * r_cross) * z_ind;

            const double decay = std::sqrt(1.0 - rho_ar1 * rho_ar1);
            res_tmax = rho_ar1 * res_tmax + decay * z1;
            res_tmin = rho_ar1 * res_tmin + decay * z2;

            const double pw = std::max(0.01, std::min(0.99, pwd / std::max(0.01, 1.0 - pww + pwd)));
            const double delta_t = (units == 1) ? 3.6 : 2.0;
            const double t_shift = is_wet ? (-(1.0 - pw) * delta_t) : (pw * delta_t);
            double tmax = mean_tmax + t_shift + sd_tmax * res_tmax;
            double tmin = mean_tmin + sd_tmin * res_tmin;
            if (tmin > tmax) std::swap(tmin, tmax);

            // 3. Precipitation depth (Pearson Type III via Wilson-Hilferty transformation)
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

            ts_prcp[k] = (dt > 0.0) ? (prcp / dt) : prcp; // mm/d or in/d
            ts_tmax[k] = tmax;                            // deg C or deg F
            ts_tmin[k] = tmin;                            // deg C or deg F
        }

        // Package and write the 3 Time Series Definitions into outargs
        // Output 1: Precip_Rate (instantaneous)
        write_time_series_record(&outargs[0], 1, 0, num_points, ts_times.data(), ts_prcp.data());
        // Output 2: Tmax (instantaneous)
        write_time_series_record(&outargs[single_ts_size], 1, 0, num_points, ts_times.data(), ts_tmax.data());
        // Output 3: Tmin (instantaneous)
        write_time_series_record(&outargs[2 * single_ts_size], 1, 0, num_points, ts_times.data(), ts_tmin.data());

        if (g_step_session.log_file && g_step_session.log_level >= LOG_INFO) {
            fprintf(g_step_session.log_file,
                "[CligenWeather] Realization %ld generated: %d points from OLE %.1f to %.1f (Total precip=%.1f %s, avg Tmax=%.1f, avg Tmin=%.1f)\n",
                real_id, num_points, ts_times.front(), ts_times.back(),
                std::accumulate(ts_prcp.begin(), ts_prcp.end(), 0.0), (units == 1 ? "in" : "mm"),
                std::accumulate(ts_tmax.begin(), ts_tmax.end(), 0.0) / num_points,
                std::accumulate(ts_tmin.begin(), ts_tmin.end(), 0.0) / num_points);
            fflush(g_step_session.log_file);
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

// ===========================================================================
// Exported GoldSim Entry Points
// ===========================================================================

// Primary entry point for GoldSim: Generates complete Time Series once per realization
GS_EXPORT void GS_CALL CligenWeather(int methodID, int* status,
                                     double* inargs, double* outargs) {
    execute_cligen_timeseries(methodID, status, inargs, outargs);
}

// Explicit alias for Time Series mode
GS_EXPORT void GS_CALL cligen_weather_ts(int methodID, int* status,
                                        double* inargs, double* outargs) {
    execute_cligen_timeseries(methodID, status, inargs, outargs);
}

// Step-by-step array mode (112 inputs: 4 scalars + 9 vectors of 12)
GS_EXPORT void GS_CALL cligen_weather_arrays(int methodID, int* status,
                                            double* inargs, double* outargs) {
    execute_cligen(true, methodID, status, inargs, outargs);
}

// Step-by-step scalar mode (13 scalar inputs)
GS_EXPORT void GS_CALL cligen_weather_step(int methodID, int* status,
                                           double* inargs, double* outargs) {
    execute_cligen(false, methodID, status, inargs, outargs);
}
