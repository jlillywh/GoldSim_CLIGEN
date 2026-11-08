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
        if (s == "OFF" || s == "NONE" || s == "0") return LOG_OFF;
        if (s == "DEBUG" || s == "3") return LOG_DEBUG;
        if (s == "INFO"  || s == "2") return LOG_INFO;
        if (s == "ERROR" || s == "1") return LOG_ERROR;
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

#include <cstdarg>

static void log_cligen_level(LogLevel req_level, const char* fmt, ...);

#define log_cligen(...)      log_cligen_level(LOG_DEBUG, __VA_ARGS__)
#define log_cligen_err(...)  log_cligen_level(LOG_ERROR, __VA_ARGS__)
#define log_cligen_info(...) log_cligen_level(LOG_INFO,  __VA_ARGS__)

namespace {
    struct StepSession {
        bool initialized = false;
        long current_realization = -1;
        long step_count = 0;
        bool prev_wet = false;      // wet/dry Markov state carried from yesterday
        double res_tmax = 0.0;      // AR(1) max-temperature residual carried forward
        double res_tmin = 0.0;      // AR(1) min-temperature residual carried forward
        long active_day_number = -999999999L;
        cligen::DayWeather cur_day_weather{};
        cligen::SubdailyStorm cur_day_storm{};
        double prev_day_tmax = 20.0;
        double prev_day_prcp = 0.0;
        cligen::SubdailyStorm prev_day_storm{};
        bool has_prev_day = false;
        cligen::MonthlyCubicSpline param_spline[NUM_PARAM_BLOCKS];
        std::mt19937_64 rng;
        std::mt19937_64 subdaily_rng;
        std::normal_distribution<double> norm{0.0, 1.0};
        std::uniform_real_distribution<double> uni{0.0, 1.0};
        LogLevel log_level = LOG_ERROR;

        void open_log() {
            log_level = detect_log_level();
            if (log_level < LOG_INFO) return;
            log_cligen_info("====================================================================================\n");
            log_cligen_info(" CligenWeather.dll v%.1f Log (Level: %s)\n", DLL_VERSION,
                log_level == LOG_DEBUG ? "DEBUG" : (log_level == LOG_INFO ? "INFO" : "ERROR"));
            log_cligen_info("====================================================================================\n");
        }

        void log_error(const char* msg) {
            log_cligen_err("[ERROR] Step %ld: %s\n", step_count, msg);
        }

        void close_log() {
            if (step_count > 0) {
                log_cligen_info("[CligenWeather] Simulation finished. Total steps: %ld\n", step_count);
            } else {
                log_cligen_info("[CligenWeather] Simulation cleanup complete.\n");
            }
        }

        void reset_realization(long real_id) {
            current_realization = real_id;
            uint64_t seed = (static_cast<uint64_t>(real_id) * 6364136223846793005ULL + 1442695040888963407ULL) ^ 0x9E3779B97F4A7C15ULL;
            if (seed == 0) seed = 189;
            rng.seed(seed);
            uint64_t sub_seed = (seed ^ 0xD1B54A32D192ED03ULL) * 0x94D049BB133111EBULL;
            if (sub_seed == 0) sub_seed = 7919;
            subdaily_rng.seed(sub_seed);
            norm.reset();
            uni.reset();
            prev_wet = false;
            res_tmax = 0.0;
            res_tmin = 0.0;
            active_day_number = -999999999L;
            prev_day_prcp = 0.0;
            prev_day_storm = cligen::SubdailyStorm{};
            has_prev_day = false;
            log_cligen_info("\n--- Realization %ld (Seed: %llu) ---\n", real_id, (unsigned long long)seed);
        }
    };

    StepSession g_step_session;
}

static void log_cligen_level(LogLevel req_level, const char* fmt, ...) {
    if (g_step_session.log_level < req_level) return;
    char buf[4096];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);

    // 1. C:\temp\cligen_debug.log
    FILE* f1 = fopen("C:\\temp\\cligen_debug.log", "a");
    if (f1) { fputs(buf, f1); fflush(f1); fclose(f1); }

    // 2. DLL folder cligen_debug.log
    std::string dll_dir = get_dll_directory();
    if (!dll_dir.empty()) {
        std::string p2 = dll_dir + "\\cligen_debug.log";
        FILE* f2 = fopen(p2.c_str(), "a");
        if (f2) { fputs(buf, f2); fflush(f2); fclose(f2); }
    }
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
        if (g_step_session.log_level >= LOG_INFO) {
            log_cligen("[CligenWeather] Session initialized (%s mode).\n",
                array_mode ? "112-vector spline" : "13-scalar step");
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
        const double raw_dt    = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        const double dt        = (raw_dt > 100.0) ? (raw_dt / 86400.0) : raw_dt;
        const long   real_id   = (inargs[2] >= 1.0) ? static_cast<long>(inargs[2]) : 1;
        const int    units     = (inargs[3] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        // Carry state forward between days; reseed only when a NEW realization begins
        if (!g_step_session.initialized || g_step_session.current_realization != real_id) {
            g_step_session.open_log();
            g_step_session.initialized = true;
            g_step_session.reset_realization(real_id);
        }

        const bool is_subdaily = (dt < 0.999);
        const double dt_hours  = std::max(1.0 / 3600.0, dt * 24.0);
        const long day_number  = static_cast<long>(std::floor(time_val + 1e-9));
        const double frac_day  = std::max(0.0, std::min(1.0, time_val - static_cast<double>(day_number)));
        const double h0        = frac_day * 24.0;
        const double h1        = std::min(24.0, h0 + dt_hours);

        int month_idx = 0;
        double pww = 0.45, pwd = 0.20, prcp_mean = 6.0, prcp_sd = 7.0, prcp_skew = 2.0;
        double mean_tmax = 20.0, sd_tmax = 4.0, mean_tmin = 10.0, sd_tmin = 3.5;

        if (!is_subdaily || day_number != g_step_session.active_day_number) {
            if (array_mode) {
                int doy = 1;
                bool leap = false;
                if (time_val >= 10000.0) {
                    int year = 0, month = 0, day = 0;
                    cligen::ole_to_date(static_cast<double>(day_number), year, month, day, doy);
                    leap = cligen::is_leap_year(year);
                    month_idx = std::max(0, std::min(11, month - 1));
                } else {
                    long day_idx = std::max(0L, day_number);
                    doy = static_cast<int>((day_idx % 365) + 1);
                    month_idx = cligen::doy_to_month(doy, false);
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

            cligen::DayClimate climate;
            climate.pww = pww;
            climate.pwd = pwd;
            climate.prcp_mean = prcp_mean;
            climate.prcp_sd = prcp_sd;
            climate.prcp_skew = prcp_skew;
            climate.mean_tmax = mean_tmax;
            climate.sd_tmax = sd_tmax;
            climate.mean_tmin = mean_tmin;
            climate.sd_tmin = sd_tmin;
            climate.units = units;

            if (g_step_session.has_prev_day) {
                g_step_session.prev_day_tmax = g_step_session.cur_day_weather.tmax;
                g_step_session.prev_day_prcp = g_step_session.cur_day_weather.prcp;
                g_step_session.prev_day_storm = g_step_session.cur_day_storm;
            } else if (is_subdaily) {
                // Draw pre-start day (d = -1) on subdaily_rng for cross-midnight spillover at t = 0
                cligen::DayMemory pre_mem{};
                std::normal_distribution<double> pre_norm(0.0, 1.0);
                cligen::DayWeather pre_day = cligen::draw_day(
                    climate, pre_mem, g_step_session.subdaily_rng, pre_norm, g_step_session.uni);
                g_step_session.prev_day_tmax = pre_day.tmax;
                g_step_session.prev_day_prcp = pre_day.prcp;
                if (pre_day.is_wet && pre_day.prcp > 0.0) {
                    g_step_session.prev_day_storm = cligen::draw_cligen_storm(
                        pre_day.prcp, units, month_idx, dt_hours,
                        g_step_session.subdaily_rng, g_step_session.uni);
                } else {
                    g_step_session.prev_day_storm = cligen::SubdailyStorm{};
                }
            }

            cligen::DayMemory mem;
            mem.prev_wet = g_step_session.prev_wet;
            mem.res_tmax = g_step_session.res_tmax;
            mem.res_tmin = g_step_session.res_tmin;
            g_step_session.cur_day_weather = cligen::draw_day(
                climate, mem, g_step_session.rng, g_step_session.norm, g_step_session.uni);
            g_step_session.prev_wet = mem.prev_wet;
            g_step_session.res_tmax = mem.res_tmax;
            g_step_session.res_tmin = mem.res_tmin;

            if (!g_step_session.has_prev_day) {
                if (!is_subdaily) {
                    g_step_session.prev_day_tmax = g_step_session.cur_day_weather.tmax;
                }
                g_step_session.has_prev_day = true;
            }
            g_step_session.active_day_number = day_number;

            if (is_subdaily && g_step_session.cur_day_weather.is_wet && g_step_session.cur_day_weather.prcp > 0.0) {
                g_step_session.cur_day_storm = cligen::draw_cligen_storm(
                    g_step_session.cur_day_weather.prcp, units, month_idx, dt_hours,
                    g_step_session.subdaily_rng, g_step_session.uni);
            } else {
                g_step_session.cur_day_storm = cligen::SubdailyStorm{};
            }
        }

        const cligen::DayWeather& day = g_step_session.cur_day_weather;
        const bool is_wet = day.is_wet;

        if (!is_subdaily) {
            outargs[0] = (dt > 0.0) ? (day.prcp / dt) : day.prcp; // Precip_Rate (mm/day or in/day)
            outargs[1] = day.tmax;                                // Tmax
            outargs[2] = day.tmin;                                // Tmin
        } else {
            const double step_depth =
                cligen::eval_storm_interval_depth(day.prcp, g_step_session.cur_day_storm, h0, h1) +
                cligen::eval_storm_interval_depth(g_step_session.prev_day_prcp, g_step_session.prev_day_storm, h0 + 24.0, h1 + 24.0);
            const double step_rate = step_depth / dt;
            const double t_diurnal = cligen::eval_diurnal_temp(
                h0, g_step_session.prev_day_tmax, day.tmin, day.tmax, day.tmin);
            outargs[0] = step_rate;
            outargs[1] = t_diurnal;
            outargs[2] = t_diurnal;
        }

        g_step_session.step_count++;
        if (g_step_session.log_level >= LOG_DEBUG) {
            if (g_step_session.step_count <= 100 || (g_step_session.step_count % 365 == 0)) {
                log_cligen(
                    "Step %5ld | Time=%6.1f | Real=%ld | In: PWW=%.3f PWD=%.3f MP=%.3f SDP=%.3f Skew=%.2f MTmax=%.1f SDTmax=%.2f MTmin=%.1f SDTmin=%.2f | Out: Wet=%d Precip_Rate=%7.4f %s Tmax=%5.2f Tmin=%5.2f\n",
                    g_step_session.step_count, time_val, real_id,
                    pww, pwd, prcp_mean, prcp_sd, prcp_skew,
                    mean_tmax, sd_tmax, mean_tmin, sd_tmin,
                    (is_wet ? 1 : 0), outargs[0], (units == 1 ? "in/d" : "mm/d"),
                    outargs[1], outargs[2]);
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

    // Log all incoming method calls to both C:\temp\cligen_debug.log and DLL folder cligen_debug.log
    log_cligen("[GoldSim Call] MethodID=%d (%s), status=%d, inargs=%p, outargs=%p\n",
               methodID,
               methodID == 0 ? "XF_INITIALIZE" :
               methodID == 1 ? "XF_CALCULATE" :
               methodID == 2 ? "XF_REP_VERSION" :
               methodID == 3 ? "XF_REP_ARGUMENTS" :
               methodID == 99 ? "XF_CLEANUP" : "UNKNOWN",
               status ? *status : -999,
               (void*)inargs, (void*)outargs);

    static long s_cached_real_id = -1;
    static double s_cached_dt_days = -1.0;
    static double s_cached_sim_dur = -1.0;
    static std::vector<double> s_cached_ts_output;

    switch (methodID) {
    case XF_INITIALIZE: {
        s_cached_real_id = -1;
        s_cached_dt_days = -1.0;
        s_cached_sim_dur = -1.0;
        s_cached_ts_output.clear();
        g_step_session.open_log();
        g_step_session.initialized = true;
        g_step_session.current_realization = -1;
        g_step_session.step_count = 0;
        log_cligen("[CligenWeather_TS] Time Series Generator initialized.\n");
        break;
    }

    case XF_REP_VERSION: {
        if (outargs) outargs[0] = DLL_VERSION;
        log_cligen("  XF_REP_VERSION: reported version %.1f\n", DLL_VERSION);
        break;
    }

    case XF_REP_ARGUMENTS: {
        if (outargs) {
            outargs[0] = static_cast<double>(NUM_TS_INPUTS); // 113 inputs
            outargs[1] = static_cast<double>(g_ts_output_buffer_size); // 500,000 doubles initial buffer
        }
        log_cligen("  XF_REP_ARGUMENTS: reported inputs=%.0f, output_buffer=%.0f\n",
                   outargs ? outargs[0] : 0.0, outargs ? outargs[1] : 0.0);
        break;
    }

    case XF_CALCULATE: {
        if (!inargs || !outargs) {
            g_step_session.log_error("Null inargs or outargs pointer passed by GoldSim");
            if (status) *status = XF_FAILURE;
            return;
        }

        const double raw_start_time = inargs[0]; // GoldSim StartTime / DateTime run property
        const double raw_sim_dur    = (inargs[1] > 0.0) ? inargs[1] : 1.0;
        const double raw_dt         = (inargs[2] > 0.0) ? inargs[2] : 1.0;
        const long   real_id        = (inargs[3] >= 1.0) ? static_cast<long>(std::round(inargs[3])) : 1;
        const int    units          = (inargs[4] >= 0.5) ? 1 : 0; // 0 = Metric, 1 = US Customary

        // Auto-detect and normalize units for SimDuration and TimestepLength:
        const bool dur_in_seconds = (raw_sim_dur > 10000.0);
        double sim_dur_days = dur_in_seconds ? (raw_sim_dur / 86400.0) : raw_sim_dur;
        double dt_days      = (raw_dt > 100.0 || (dur_in_seconds && raw_dt > 1.0)) ? (raw_dt / 86400.0) : raw_dt;
        if (dt_days <= 0.0) dt_days = 1.0;
        if (sim_dur_days <= 0.0) sim_dur_days = 1.0;

        int num_steps = static_cast<int>(std::ceil(sim_dur_days / dt_days));
        if (num_steps < 1) num_steps = 1;

        // GoldSim requires: "The last data point must be on or after the simulation end time."
        constexpr int BUFFER_STEPS = 5;
        int num_points = num_steps + 1 + BUFFER_STEPS;

        int single_ts_size = 8 + 2 * num_points;
        int total_req_size = NUM_TS_OUTPUTS * single_ts_size;

        if (total_req_size > g_ts_output_buffer_size) {
            outargs[0] = static_cast<double>(total_req_size);
            g_ts_output_buffer_size = total_req_size;
            if (status) *status = -2; // XF_INCREASE_MEMORY
            return;
        }

        // Realization skip guard:
        // If GoldSim triggers XF_CALCULATE again during the same realization with identical settings,
        // instantly return the cached time series without re-running stochastic generation.
        if (s_cached_real_id == real_id &&
            std::fabs(s_cached_dt_days - dt_days) < 1e-12 &&
            std::fabs(s_cached_sim_dur - sim_dur_days) < 1e-9 &&
            s_cached_ts_output.size() == static_cast<size_t>(total_req_size)) {
            std::memcpy(outargs, s_cached_ts_output.data(), total_req_size * sizeof(double));
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

        // Auto-detect Calendar vs Elapsed mode:
        // In GoldSim Time Series Definitions:
        // - Calendar mode (time_basis=1.0): GoldSim requires timestamps in SECONDS (Julian Day * 86,400.0).
        // - Elapsed mode (time_basis=0.0): GoldSim requires timestamps in elapsed days (k * dt_days).
        const bool is_calendar = (raw_start_time >= 10000.0);
        const double time_basis = is_calendar ? 1.0 : 0.0;
        constexpr double SEC_PER_DAY = 86400.0;

        const bool is_subdaily = (dt_days < 0.999);
        const double dt_hours  = std::max(1.0 / 3600.0, dt_days * 24.0);
        const double start_day_floor = is_calendar ? std::floor(raw_start_time + 1e-12) : 0.0;
        const double initial_frac_day = is_calendar ? std::max(0.0, raw_start_time - start_day_floor) : 0.0;
        const int max_day_idx = std::max(0, static_cast<int>(std::floor(initial_frac_day + (num_points - 1) * dt_days + 1e-9)));
        const int num_sim_days = max_day_idx + 2; // +1 lookahead day for continuous midnight diurnal cooling

        std::vector<cligen::DayWeather> sim_days(num_sim_days);
        std::vector<cligen::SubdailyStorm> sim_storms(num_sim_days);
        cligen::DayWeather pre_day{};
        cligen::SubdailyStorm pre_storm{};

        // 1. Generate daily CLIGEN weather strictly ONCE per calendar day (daily RNG stream)
        //    and draw sub-daily storm parameters on the isolated subdaily_rng stream.
        for (int d = 0; d < num_sim_days; ++d) {
            int doy = 1;
            int month_idx = 0;
            bool leap = false;
            if (is_calendar) {
                int year = 0, month = 0, day_of_m = 0;
                cligen::ole_to_date(start_day_floor + static_cast<double>(d), year, month, day_of_m, doy);
                leap = cligen::is_leap_year(year);
                month_idx = std::max(0, std::min(11, month - 1));
            } else {
                doy = (d % 365) + 1;
                leap = false;
                month_idx = cligen::doy_to_month(doy, false);
            }

            const double u = cligen::month_phase(doy, leap);
            double param[NUM_PARAM_BLOCKS];
            for (int blk = 0; blk < NUM_PARAM_BLOCKS; ++blk) {
                param[blk] = param_spline[blk].evaluate(u);
            }

            cligen::DayClimate climate;
            climate.pww = std::max(0.01, std::min(0.99, param[BLK_PWW]));
            climate.pwd = std::max(0.01, std::min(0.99, param[BLK_PWD]));
            climate.prcp_mean = param[BLK_MEANP];
            climate.prcp_sd = param[BLK_SDP];
            climate.prcp_skew = std::max(0.0, std::min(4.0, param[BLK_SKEWP]));
            climate.mean_tmax = param[BLK_MTMAX];
            climate.sd_tmax = param[BLK_SDTMAX];
            climate.mean_tmin = param[BLK_MTMIN];
            climate.sd_tmin = param[BLK_SDTMIN];
            climate.units = units;

            if (d == 0 && is_subdaily) {
                // Draw pre-start day (d = -1) on subdaily_rng so storms starting late on Day -1
                // naturally spill across the initial midnight (t = 0) without perturbing daily rng.
                cligen::DayMemory pre_mem{};
                std::normal_distribution<double> pre_norm(0.0, 1.0);
                pre_day = cligen::draw_day(
                    climate, pre_mem, g_step_session.subdaily_rng, pre_norm, g_step_session.uni);
                if (pre_day.is_wet && pre_day.prcp > 0.0) {
                    pre_storm = cligen::draw_cligen_storm(
                        pre_day.prcp, units, month_idx, dt_hours,
                        g_step_session.subdaily_rng, g_step_session.uni);
                }
            }

            cligen::DayMemory mem;
            mem.prev_wet = prev_wet;
            mem.res_tmax = res_tmax;
            mem.res_tmin = res_tmin;
            sim_days[d] = cligen::draw_day(
                climate, mem, g_step_session.rng, g_step_session.norm, g_step_session.uni);
            prev_wet = mem.prev_wet;
            res_tmax = mem.res_tmax;
            res_tmin = mem.res_tmin;

            if (is_subdaily && sim_days[d].is_wet && sim_days[d].prcp > 0.0) {
                sim_storms[d] = cligen::draw_cligen_storm(
                    sim_days[d].prcp, units, month_idx, dt_hours,
                    g_step_session.subdaily_rng, g_step_session.uni);
            }
        }

        // 2. Populate Time Series output points (daily values when dt >= 1 day;
        //    CLIGEN double-exponential storm + diurnal temperature on both Tmax & Tmin when dt < 1 day).
        for (int k = 0; k < num_points; ++k) {
            const double cur_day = raw_start_time + k * dt_days;
            ts_times[k] = is_calendar ? (cur_day * SEC_PER_DAY) : (k * dt_days);

            const double rel_days = initial_frac_day + k * dt_days;
            const int d = std::max(0, std::min(num_sim_days - 2, static_cast<int>(std::floor(rel_days + 1e-9))));
            const cligen::DayWeather& day = sim_days[d];

            if (!is_subdaily) {
                ts_prcp[k] = (dt_days > 0.0) ? (day.prcp / dt_days) : day.prcp; // mm/d or in/d
                ts_tmax[k] = day.tmax;                                          // deg C or deg F
                ts_tmin[k] = day.tmin;                                          // deg C or deg F
            } else {
                const double h0 = std::max(0.0, std::min(24.0, (rel_days - static_cast<double>(d)) * 24.0));
                const double h1 = std::min(24.0, h0 + dt_hours);

                const double prev_prcp = (d > 0) ? sim_days[d - 1].prcp : pre_day.prcp;
                const cligen::SubdailyStorm& prev_st = (d > 0) ? sim_storms[d - 1] : pre_storm;
                const double step_depth =
                    cligen::eval_storm_interval_depth(day.prcp, sim_storms[d], h0, h1) +
                    cligen::eval_storm_interval_depth(prev_prcp, prev_st, h0 + 24.0, h1 + 24.0);
                const double step_rate = step_depth / dt_days;

                const double prev_tmax = (d > 0) ? sim_days[d - 1].tmax : day.tmax;
                const double next_tmin = sim_days[d + 1].tmin;
                const double t_diurnal = cligen::eval_diurnal_temp(h0, prev_tmax, day.tmin, day.tmax, next_tmin);

                ts_prcp[k] = step_rate;
                ts_tmax[k] = t_diurnal;
                ts_tmin[k] = t_diurnal;
            }
        }

        // Package and write the 3 Time Series Definitions into outargs
        // Output 1: Precip_Rate (instantaneous)
        write_time_series_record(&outargs[0], static_cast<int>(time_basis), 0, num_points, ts_times.data(), ts_prcp.data());
        // Output 2: Tmax (instantaneous)
        write_time_series_record(&outargs[single_ts_size], static_cast<int>(time_basis), 0, num_points, ts_times.data(), ts_tmax.data());
        // Output 3: Tmin (instantaneous)
        write_time_series_record(&outargs[2 * single_ts_size], static_cast<int>(time_basis), 0, num_points, ts_times.data(), ts_tmin.data());

        // Detailed diagnostic logging to both C:\temp\cligen_debug.log and local cligen_debug.log
        log_cligen("=======================================================================\n");
        log_cligen(" CLIGEN Time Series Generation — Realization %ld\n", real_id);
        log_cligen("=======================================================================\n");
        log_cligen("Raw inargs:  StartTime=%.6f, SimDur=%.6f, TimeStep=%.6f, Real=%ld, Units=%d\n",
                   raw_start_time, raw_sim_dur, raw_dt, real_id, units);
        log_cligen("Normalized:  SimDurDays=%.2f, DtDays=%.4f, NumSteps=%d, NumPoints=%d (+%d buffer)\n",
                   sim_dur_days, dt_days, num_steps, num_points, BUFFER_STEPS);
        log_cligen("Time Basis:  %s (code=%.0f), t_first=%.1f s (day %.4f), t_last=%.1f s (day %.4f)\n",
                   is_calendar ? "Calendar" : "Elapsed", time_basis,
                   ts_times.front(), is_calendar ? (ts_times.front() / SEC_PER_DAY) : ts_times.front(),
                   ts_times.back(), is_calendar ? (ts_times.back() / SEC_PER_DAY) : ts_times.back());
        if (is_calendar) {
            int y1, m1, d1, doy1, y2, m2, d2, doy2;
            cligen::ole_to_date(ts_times.front() / SEC_PER_DAY, y1, m1, d1, doy1);
            cligen::ole_to_date(ts_times.back() / SEC_PER_DAY, y2, m2, d2, doy2);
            log_cligen("Calendar:    Start=%04d-%02d-%02d to End=%04d-%02d-%02d\n",
                       y1, m1, d1, y2, m2, d2);
        }
        log_cligen("Output 1 (Precip_Rate) @ outargs[0]:\n");
        log_cligen("  Header: Format=%.0f, Ver=%.0f, Basis=%.0f, Type=%.0f, Rows=%.0f, Cols=%.0f, Series=%.0f, Pts=%.0f\n",
                   outargs[0], outargs[1], outargs[2], outargs[3], outargs[4], outargs[5], outargs[6], outargs[7]);
        log_cligen("  First 3 points: (t=%.1f, v=%.3f), (t=%.1f, v=%.3f), (t=%.1f, v=%.3f)\n",
                   outargs[8], outargs[8 + num_points],
                   outargs[9], outargs[8 + num_points + 1],
                   outargs[10], outargs[8 + num_points + 2]);
        log_cligen("  Last 3 points:  (t=%.1f, v=%.3f), (t=%.1f, v=%.3f), (t=%.1f, v=%.3f)\n",
                   outargs[8 + num_points - 3], outargs[8 + 2 * num_points - 3],
                   outargs[8 + num_points - 2], outargs[8 + 2 * num_points - 2],
                   outargs[8 + num_points - 1], outargs[8 + 2 * num_points - 1]);

        log_cligen("Output 2 (Tmax) @ outargs[%d]:\n", single_ts_size);
        log_cligen("  Header: Format=%.0f, Ver=%.0f, Basis=%.0f, Type=%.0f, Rows=%.0f, Cols=%.0f, Series=%.0f, Pts=%.0f\n",
                   outargs[single_ts_size + 0], outargs[single_ts_size + 1], outargs[single_ts_size + 2],
                   outargs[single_ts_size + 3], outargs[single_ts_size + 4], outargs[single_ts_size + 5],
                   outargs[single_ts_size + 6], outargs[single_ts_size + 7]);

        log_cligen("Output 3 (Tmin) @ outargs[%d]:\n", 2 * single_ts_size);
        log_cligen("  Header: Format=%.0f, Ver=%.0f, Basis=%.0f, Type=%.0f, Rows=%.0f, Cols=%.0f, Series=%.0f, Pts=%.0f\n",
                   outargs[2 * single_ts_size + 0], outargs[2 * single_ts_size + 1], outargs[2 * single_ts_size + 2],
                   outargs[2 * single_ts_size + 3], outargs[2 * single_ts_size + 4], outargs[2 * single_ts_size + 5],
                   outargs[2 * single_ts_size + 6], outargs[2 * single_ts_size + 7]);
        log_cligen("Total doubles written: %d (offset 0 to %d)\n", total_req_size, total_req_size - 1);

        // Store generated series in memory cache
        s_cached_ts_output.assign(outargs, outargs + total_req_size);
        s_cached_real_id = real_id;
        s_cached_dt_days = dt_days;
        s_cached_sim_dur = sim_dur_days;
        break;
    }

    case XF_CLEANUP: {
        s_cached_real_id = -1;
        s_cached_dt_days = -1.0;
        s_cached_sim_dur = -1.0;
        s_cached_ts_output.clear();
        g_step_session.close_log();
        g_step_session.initialized = false;
        log_cligen("  XF_CLEANUP: closed.\n");
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

// Primary entry point: 113 inputs (5 scalars + 9 monthly 12-item vectors) -> 3 Time Series Definition outputs
GS_EXPORT void GS_CALL CligenWeather(int methodID, int* status,
                                     double* inargs, double* outargs) {
    execute_cligen_timeseries(methodID, status, inargs, outargs);
}

// Explicit alias for GoldSim Time Series Definition mode (113 inputs -> 3 Time Series)
GS_EXPORT void GS_CALL cligen_weather_ts(int methodID, int* status,
                                        double* inargs, double* outargs) {
    execute_cligen_timeseries(methodID, status, inargs, outargs);
}

// Step-by-step mode alias: 112 inputs (4 scalars + 9 monthly 12-item vectors) -> 3 daily scalar outputs
GS_EXPORT void GS_CALL cligen_weather_step(int methodID, int* status,
                                           double* inargs, double* outargs) {
    execute_cligen(true, methodID, status, inargs, outargs);
}

// Explicit alias for 112-input monthly vector spline mode
GS_EXPORT void GS_CALL cligen_weather_arrays(int methodID, int* status,
                                            double* inargs, double* outargs) {
    execute_cligen(true, methodID, status, inargs, outargs);
}

// Legacy daily scalar mode (13 scalar inputs)
GS_EXPORT void GS_CALL cligen_weather_scalar(int methodID, int* status,
                                             double* inargs, double* outargs) {
    execute_cligen(false, methodID, status, inargs, outargs);
}
