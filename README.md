# GoldSim - CLIGEN Integration

[![Windows 64-bit](https://img.shields.io/badge/platform-Windows%2064--bit-blue.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[GoldSim](https://www.goldsim.com/Web/Home/) DLL `CligenWeather.dll` links directly to GoldSim's External Element to run stochastic simulations of daily weather data (precipitation, Tmax, and Tmin).

---

## Why CLIGEN vs. WGEN?

Both **CLIGEN** (USDA-ARS) and **WGEN** (Richardson & Wright) are established stochastic weather generators that share a first-order two-state Markov chain for wet/dry day transitions. However, CLIGEN was developed specifically to overcome three key limitations in WGEN:

1. **Extreme Precipitation Tails (Pearson Type III vs. Gamma)**:
   * **WGEN** samples daily precipitation depth from a 2-parameter Gamma or Exponential distribution. Because these distributions have exponentially thin tails and no skewness adjustment, WGEN systematically **underpredicts extreme storm depths** (25-year to 100-year events) by 20% to 50%.
   * **CLIGEN** uses a 3-parameter **Pearson Type III distribution** via the Wilson-Hilferty transformation. By incorporating monthly skewness coefficients ($\text{Skew}_P$), CLIGEN accurately reproduces both everyday moderate rainfall and the heavy positive tails of extreme design storm events.

2. **Storm Patterning & Sub-Daily Mechanics**:
   * **WGEN** produces only bulk daily totals without storm pattern parameterization.
   * **CLIGEN** incorporates internal sub-daily mechanics developed for the USDA Water Erosion Prediction Project (WEPP), including time-to-peak and peak intensity modeling using a mass-conserving double-exponential profile, and diurnal temperature cycles conditioned on daily $T_{\max}$ and $T_{\min}$.

3. **Physically Conditioned Temperatures**:
   * **WGEN** uses an unconditioned continuous autoregressive model that can occasionally produce runaway temperature extremes (e.g., $>55^\circ\text{C}$ / $130^\circ\text{F}$) over long Monte Carlo realizations.
   * **CLIGEN** conditions temperature directly on wet/dry occurrence: on wet days, $T_{\max}$ is depressed to reflect cloud cover and evaporative cooling. Coupled with bivariate cross-correlation ($r = 0.65$) and lag-1 persistence ($\rho = 0.60$), it prevents runaway extremes while preserving day-to-day temperature memory.

---

## Technical Interface Specification

The DLL exports two entry points to accommodate either **smooth seasonal cubic spline interpolation** (passing 12-month arrays) or **direct daily scalar values**:

### 1. Vector Spline Mode: `CligenWeather` or `cligen_weather_arrays` (Recommended)
Passes the full 12-month climate tables to the DLL. The DLL fits a **periodic cubic spline** through each 12-month parameter cycle, eliminating artificial step-function jumps at month boundaries and maintaining continuous first and second derivatives across the year.

* **Function Name**: `CligenWeather` (or `cligen_weather_arrays`)
* **Calling Convention**: `__cdecl`
* **Total GoldSim Input Arguments**: **112** (4 scalar controls + 9 vectors of 12 elements each)

| Index Range | Argument Name | GoldSim Type | Description |
| :---: | :--- | :--- | :--- |
| `[0]` | `Simulation_Time` | Scalar Value (`day`) | Elapsed simulation time (GoldSim `ETIME`) |
| `[1]` | `Time_Step` | Scalar Value (`day`) | Timestep duration (GoldSim `DT`, e.g. 1.0 day) |
| `[2]` | `Realization` | Scalar Value (`-`) | Current Monte Carlo realization index |
| `[3]` | `Units_Flag` | Scalar Value (`-`) | `0` = Metric (mm, &deg;C), `1` = US Customary (in, &deg;F) |
| `[4..15]` | `PWW[12]` | 12-item Vector | Monthly $P(W\|W)$ wet-given-wet probabilities (Jan..Dec) |
| `[16..27]` | `PWD[12]` | 12-item Vector | Monthly $P(W\|D)$ wet-given-dry probabilities (Jan..Dec) |
| `[28..39]` | `Mean_P[12]` | 12-item Vector | Monthly wet-day mean precipitation depth (`mm` or `in`) |
| `[40..51]` | `SD_P[12]` | 12-item Vector | Monthly daily precip standard deviation (`mm` or `in`) |
| `[52..63]` | `Skew_P[12]` | 12-item Vector | Monthly daily precip skewness coefficient |
| `[64..75]` | `Mean_Tmax[12]` | 12-item Vector | Monthly mean daily max air temperature (`&deg;C` or `&deg;F`) |
| `[76..87]` | `SD_Tmax[12]` | 12-item Vector | Monthly max temp daily standard deviation (`Cdeg` or `Fdeg`) |
| `[88..99]` | `Mean_Tmin[12]` | 12-item Vector | Monthly mean daily min air temperature (`&deg;C` or `&deg;F`) |
| `[100..111]` | `SD_Tmin[12]` | 12-item Vector | Monthly min temp daily standard deviation (`Cdeg` or `Fdeg`) |

> **Note on GoldSim Arguments:** When linking 9 vectors of dimension 12 alongside the 4 control scalars, GoldSim reports $4 + 9 \times 12 = 112$ scalar arguments to the external function. The DLL reconstructs day-of-year from `ETIME` and `DT` and evaluates each spline smoothly at mid-month coordinates.

---

### 2. Daily Scalar Mode: `cligen_weather_step`
For models where you prefer to compute or lookup single daily scalar parameters directly inside GoldSim and pass only 13 scalars:

* **Function Name**: `cligen_weather_step`
* **Calling Convention**: `__cdecl`
* **Total GoldSim Input Arguments**: **13** Scalars (`Simulation_Time`, `Time_Step`, `Realization`, `Units_Flag`, `P_WW`, `P_WD`, `Mean_P`, `SD_P`, `Skew_P`, `Mean_Tmax`, `SD_Tmax`, `Mean_Tmin`, `SD_Tmin`).

---

### Outputs (`outargs` &mdash; 3 Scalars)
Both entry points produce the same 3 daily output variables:
| Index | Output Parameter | Units | Description |
| :---: | :--- | :--- | :--- |
| `[0]` | `Precip_Rate` | `in/day` or `mm/day` | Stochastic precipitation rate for the current timestep (`prcp / DT`) |
| `[1]` | `Tmax` | `&deg;F` or `&deg;C` | Daily maximum air temperature |
| `[2]` | `Tmin` | `&deg;F` or `&deg;C` | Daily minimum air temperature |

---

## GoldSim Setup Guide

1. **Place DLL**: Copy `CligenWeather.dll` into the folder containing your GoldSim model file (`.gsm`).
2. **Insert External Element**:
   * Right-click canvas &rarr; **Insert Element &rarr; Functions &rarr; External (DLL)**.
   * **File Name**: `CligenWeather.dll`
   * **Function**: `CligenWeather` (for 112-item monthly vector spline mode) or `cligen_weather_step` (for 13-item scalar mode).
   * **Calling Convention**: `__cdecl`
   * **Run in separate process**: Check the **"Run in separate process"** checkbox on the External Element's Definition tab. GoldSim uses 32-bit processes by default, and this option allows 32-bit GoldSim to run the 64-bit DLL via inter-process communication without Windows Error 193.
3. **Link Arguments**:
   * For **Vector Spline Mode**: Link `Simulation_Time`, `Time_Step`, `Realization`, `Units_Flag`, followed by the 9 monthly 12-item vectors in exact order.
   * For **Scalar Mode**: Link the 4 control scalars and the 9 daily scalar values.
4. **Run Simulation**:
   Run Monte Carlo realizations. The DLL automatically detects new realizations and seeds each sequence independently.

---

## Logging & Performance Control

By default, the DLL operates in **`ERROR`** mode for maximum execution speed, with zero disk I/O on timesteps.

If you wish to control or enable debug logging:
* **Default (`ERROR` / Fast Mode)**: Only critical failures (e.g. null arguments) are recorded.
* **To Enable Debug Logging (`DEBUG`)**:
  * Set the Windows environment variable: `CLIGEN_LOG_LEVEL=DEBUG`
  * **Or** simply create an empty file named `cligen_debug.txt` in the same directory as your `.gsm` model.
  * When active, detailed step progression is written to `cligen_debug.log`.
* **To Fully Disable All Logging (`OFF`)**:
  * Set environment variable: `CLIGEN_LOG_LEVEL=OFF` or `CLIGEN_LOG_LEVEL=NONE`.

---

## Building from Source

You do not need to install anything if using the pre-compiled `CligenWeather.dll`. If you wish to compile it yourself:

### Windows (MSVC or MinGW)
Double-click `build.bat` or run:
```cmd
build.bat
```

Or manually:
* **MinGW-w64 (g++)**:
  ```cmd
  g++ -O2 -std=c++17 -Wall -shared goldsim_external.cpp -o CligenWeather.dll -static -Wl,--add-stdcall-alias
  ```
* **Microsoft Visual C++ (MSVC)**:
  ```cmd
  cl /O2 /std:c++17 /EHsc /LD goldsim_external.cpp /Fe:CligenWeather.dll
  ```

### Linux (Cross-Compilation to Windows 64-bit)
```sh
make dll CXX=x86_64-w64-mingw32-g++
```

---

## License

MIT License &mdash; Free for academic, personal, and commercial engineering use.
Developed by [Lillywhite Consulting](https://lillywhiteconsulting.com).
