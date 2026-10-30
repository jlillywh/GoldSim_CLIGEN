# GoldSim – CLIGEN Native Weather Generator (`CligenWeather.dll`)

A standalone, zero-dependency **64-bit native C++17 Windows DLL (`CligenWeather.dll`)** that embeds the stochastic daily weather generation mathematics of **CLIGEN** and **WGEN** directly inside a **GoldSim External (DLL) Element**, outputting full **Time Series Definitions** at initialization (`ETime = 0`) for high-speed Monte Carlo simulations.

---

## Key Features

- **Single 64-bit Native DLL (`CligenWeather.dll`)**: Statically linked C++17 binary—no Python runtime, external libraries, or installation required.
- **Ready-to-Run Salt Lake City Example (`CLIGEN_SLC.gsm` & `CLIGEN_SLC.gsp`)**: Includes both a complete GoldSim model file (`CLIGEN_SLC.gsm`) and a GoldSim Player file (`CLIGEN_SLC.gsp`) parameterized with historical weather data from **Salt Lake City International Airport (KSLC)** so you can immediately run and inspect a working model.
- **Full Time Series Definition Output**: Computes the entire realization's daily weather record at `ETime = 0` (`XF_CALCULATE_TS`) and passes three **Time Series Definitions** (`Precip_Rate`, `Tmax`, `Tmin`) directly to GoldSim Time Series elements.
- **Seasonal Periodic Cubic Spline Interpolation**: Fits a periodic cubic spline across your 12 monthly climate parameters (`Jan..Dec`) so daily wet/dry transition probabilities, skewed precipitation depths, and temperature distributions transition smoothly across every calendar day of the year without step-function jumps at month boundaries.
- **CLIGEN Skewed Normal / Pearson Type III Precipitation & Wet/Dry Conditioned Temperatures**: Preserves high-tail extreme storm events, first-order Markov wet/dry persistence (`P_WW`, `P_WD`), and Matalas (1967) AR(1) lag-1 serial & cross-correlated daily temperature residuals (`Tmax`, `Tmin`).
- **Silent High-Speed Execution by Default**: Default log level is `ERROR` (zero disk I/O during normal runs). Create an empty file named `cligen_debug.txt` or `cligen_info.txt` in your model folder when troubleshooting is needed.

---

## Quick Start in GoldSim

### 1. Place `CligenWeather.dll` in Your Model Folder
Download `CligenWeather.dll` from the latest [GitHub Release](https://github.com/jlillywh/GoldSim_CLIGEN/releases/latest) (or from the **GoldSim DLL** tab of the Climate Simulator web app) and place it in the same folder as your `.gsm` model file.

### 2. Create 9 Monthly Vector Data Elements (`Vector[Months]`)
In GoldSim, create an Array Label Set named `Months` (12 items: `Jan..Dec`) and create **9 Data elements** defined as `1-D Array (Vector)` indexed by `Months`. Paste the 12 monthly values from the Climate Simulator **CLIGEN Parameters** tab into:

| Element Name | Units | Description |
| :--- | :--- | :--- |
| `PWW` | `-` | Probability of a wet day following a wet day (12 monthly values) |
| `PWD` | `-` | Probability of a wet day following a dry day (12 monthly values) |
| `Mean_P` | `mm` or `in` | Mean daily precipitation depth on wet days (12 monthly values) |
| `SD_P` | `mm` or `in` | Standard deviation of wet-day precipitation depth (12 monthly values) |
| `Skew_P` | `-` | Skewness coefficient of wet-day precipitation depth (12 monthly values) |
| `Mean_Tmax` | `C` or `F` | Mean daily maximum temperature (12 monthly values) |
| `SD_Tmax` | `Cdeg` or `Fdeg` | Standard deviation of daily maximum temperature (`Difference` unit type) |
| `Mean_Tmin` | `C` or `F` | Mean daily minimum temperature (12 monthly values) |
| `SD_Tmin` | `Cdeg` or `Fdeg` | Standard deviation of daily minimum temperature (`Difference` unit type) |

### 3. Configure the GoldSim External Element (`CLIGEN_DLL`)
Insert an **External** element (`Insert > Function Elements > External`) and configure its **Definition** tab:
- **DLL Path (file name)**: `CligenWeather.dll`
- **Function Name (case-sensitive)**: `CligenWeather`
- **Run in separate process**: Leave **unchecked** for 64-bit GoldSim (check only if running legacy 32-bit GoldSim).

![GoldSim External Element Interface Configuration](https://lillywhite-web-media-25197.s3.amazonaws.com/docs/captures/goldsim-external-interface-1793398167.png)

### 4. Configure the Interface Tab (14 Inputs &rarr; 3 Time Series Definition Outputs)

Add the **14 Input Arguments** (`113` double values total = 5 scalars + 9 vectors of 12) in this **exact order**:

| # | Input ID | Type & Units | GoldSim Link / Value |
| :--- | :--- | :--- | :--- |
| **1** | `Start_Time` | Scalar (`day`) | Link to GoldSim run property `StartTime` |
| **2** | `Sim_Dur` | Scalar (`day`) | Link to GoldSim run property `SimDuration` |
| **3** | `Time_Step` | Scalar (`day`) | Link to GoldSim run property `TimestepLength` (e.g., `1 day`) |
| **4** | `Realization` | Scalar (`-`) | Link to GoldSim run property `Realization` |
| **5** | `Units_Flag` | Scalar (`-`) | `0` = Metric (`mm`, `C`), `1` = US Customary (`in`, `F`) |
| **6** | `PWW` | `Vector[Months]` (`-`) | Link to `PWW` Data element (12 items) |
| **7** | `PWD` | `Vector[Months]` (`-`) | Link to `PWD` Data element (12 items) |
| **8** | `Mean_P` | `Vector[Months]` (`mm` or `in`) | Link to `Mean_P` Data element (12 items) |
| **9** | `SD_P` | `Vector[Months]` (`mm` or `in`) | Link to `SD_P` Data element (12 items) |
| **10** | `Skew_P` | `Vector[Months]` (`-`) | Link to `Skew_P` Data element (12 items) |
| **11** | `Mean_Tmax` | `Vector[Months]` (`C` or `F`) | Link to `Mean_Tmax` Data element (12 items) |
| **12** | `SD_Tmax` | `Vector[Months]` (`Cdeg` or `Fdeg`) | Link to `SD_Tmax` Data element (12 items) |
| **13** | `Mean_Tmin` | `Vector[Months]` (`C` or `F`) | Link to `Mean_Tmin` Data element (12 items) |
| **14** | `SD_Tmin` | `Vector[Months]` (`Cdeg` or `Fdeg`) | Link to `SD_Tmin` Data element (12 items) |

Add the **3 Output Arguments** as **Time Series Definitions** in this **exact order**:

| # | Output ID | Type & Units | Description |
| :--- | :--- | :--- | :--- |
| **1** | `Precip_Rate` | `Time Series Definition` (`mm/d` or `in/d`, instantaneous values) | Stochastic daily precipitation rate time series |
| **2** | `Tmax` | `Time Series Definition` (`C` or `F`, instantaneous values) | Stochastic daily maximum air temperature time series |
| **3** | `Tmin` | `Time Series Definition` (`C` or `F`, instantaneous values) | Stochastic daily minimum air temperature time series |

### 5. Link Downstream Time Series Elements
Create three **Time Series** elements in GoldSim (`Precip_DLL`, `Tmax_DLL`, `Tmin_DLL`):
1. Set **Data Source** to **Linked to external Time Series Definition**.
2. Link each Time Series element to the corresponding output from the External element:
   - `Precip_DLL` &rarr; `CLIGEN_DLL.Precip_Rate`
   - `Tmax_DLL` &rarr; `CLIGEN_DLL.Tmax`
   - `Tmin_DLL` &rarr; `CLIGEN_DLL.Tmin`

---

## Diagnostic Logging (Optional)

By default, `CligenWeather.dll` runs at `ERROR` log level and creates **no log file** unless an error occurs. To inspect initialization or daily values in `cligen_debug.log`:
- Create an empty file named `cligen_info.txt` in your GoldSim model folder for realization summary logs (`INFO` level).
- Create an empty file named `cligen_debug.txt` in your GoldSim model folder for detailed step-by-step diagnostic logs (`DEBUG` level).

---

## Building from Source

Requires `x86_64-w64-mingw32-g++` (MinGW-w64 cross-compiler):

```bash
make dll
```
