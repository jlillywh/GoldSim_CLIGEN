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

### External Element Export
* **DLL Name**: `CligenWeather.dll`
* **Function Name**: `CligenWeather` (or `cligen_weather_step` / `cligen_weather_arrays`)
* **Calling Convention**: `__cdecl` (Standard C)
* **Total GoldSim Input Arguments**: **112** (4 scalar controls + 9 monthly vectors of 12 elements each)

The DLL takes the 12-month climate parameter vectors (`Vector[12]`, Jan..Dec) and fits a **periodic cubic spline** through each 12-month parameter cycle, eliminating artificial step-function jumps at month boundaries and maintaining continuous first and second derivatives across the year.

### Inputs (`inargs` &mdash; 4 Scalars + 9 Monthly Vectors)
| Index Range | Argument Name | GoldSim Type | Typical Unit | Description |
| :---: | :--- | :--- | :--- | :--- |
| `[0]` | `Simulation_Time` | Scalar Value | `day` | GoldSim `DateTime` run property (preferred) or `ETIME` |
| `[1]` | `Time_Step` | Scalar Value | `day` | Model timestep duration (GoldSim `DT`, e.g. 1.0 day) |
| `[2]` | `Realization` | Scalar Value | `-` | Current Monte Carlo realization (GoldSim `Realization`) |
| `[3]` | `Units_Flag` | Scalar Value | `-` | `0` = Metric (mm, &deg;C), `1` = US Customary (in, &deg;F) |
| `[4..15]` | `P_WW` | Vector [12] | `-` | Monthly wet-day probability given yesterday was wet $P(W\|W)$ (12 items, Jan–Dec) |
| `[16..27]` | `P_WD` | Vector [12] | `-` | Monthly wet-day probability given yesterday was dry $P(W\|D)$ (12 items, Jan–Dec) |
| `[28..39]` | `Mean_P` | Vector [12] | `in` or `mm` | Monthly mean daily precipitation depth on wet days (12 items, Jan–Dec) |
| `[40..51]` | `SD_P` | Vector [12] | `in` or `mm` | Monthly daily precipitation standard deviation (12 items, Jan–Dec) |
| `[52..63]` | `Skew_P` | Vector [12] | `-` | Monthly daily precipitation skewness coefficient (12 items, Jan–Dec) |
| `[64..75]` | `Mean_Tmax` | Vector [12] | `&deg;F` or `&deg;C` | Monthly mean daily maximum air temperature (12 items, Jan–Dec) |
| `[76..87]` | `SD_Tmax` | Vector [12] | `Fdeg` or `Cdeg` | Monthly daily maximum temperature standard deviation (12 items, Jan–Dec) |
| `[88..99]` | `Mean_Tmin` | Vector [12] | `&deg;F` or `&deg;C` | Monthly mean daily minimum air temperature (12 items, Jan–Dec) |
| `[100..111]` | `SD_Tmin` | Vector [12] | `Fdeg` or `Cdeg` | Monthly daily minimum temperature standard deviation (12 items, Jan–Dec) |

> **Note on Calendar Synchronization:** When you pass GoldSim's built-in `DateTime` run property (e.g. `Simulation \ Run Properties \ DateTime`) as argument `[0]`, the DLL automatically detects the calendar date, extracts the day-of-year (1..366), accounts for leap years, and evaluates the splines at the correct season regardless of whether your simulation starts in January, October, or any other month. If an elapsed time `ETIME` (starting from 0) is passed, the DLL falls back to starting at Day 1 (Jan 1).

### Outputs (`outargs` &mdash; 3 Scalars)
| Index | Output Parameter | Units | Description |
| :---: | :--- | :--- | :--- |
| `[0]` | `Precip_Rate` | `in/day` or `mm/day` | Stochastic precipitation rate for the current timestep (`prcp / DT`) |
| `[1]` | `Tmax` | `&deg;F` or `&deg;C` | Daily maximum air temperature |
| `[2]` | `Tmin` | `&deg;F` or `&deg;C` | Daily minimum air temperature |

---

## GoldSim Setup Guide

1. **Place DLL**: Copy `CligenWeather.dll` into the folder containing your GoldSim model file (`.gsm`).
   * **64-bit GoldSim**: Use `CligenWeather.dll` (or `CligenWeather_x64.dll`) in-process.
   * **32-bit GoldSim**: Either use `CligenWeather_win32.dll` (renamed to `CligenWeather.dll`) in-process, or check the **"Run in separate process"** checkbox on the Definition tab when using the 64-bit `CligenWeather.dll`.
2. **Define Monthly Vector Parameters**:
   Create 12-item monthly vector Data elements (`Vector[Months]`, Jan–Dec) in GoldSim for each of the 9 climate parameters (`P_WW`, `P_WD`, `Mean_P`, `SD_P`, `Skew_P`, `Mean_Tmax`, `SD_Tmax`, `Mean_Tmin`, `SD_Tmin`).
3. **Insert External Element**:
   * Right-click canvas &rarr; **Insert Element &rarr; Functions &rarr; External (DLL)**.
   * **File Name**: `CligenWeather.dll`
   * **Function**: `CligenWeather` (or `cligen_weather_step`)
   * **Calling Convention**: `__cdecl`
4. **Link Arguments**:
   Link `Simulation_Time` (`DateTime` or `ETime`), `Time_Step` (`DT`), `Realization`, `Units_Flag`, followed by the 9 monthly 12-item vectors (`P_WW`, `P_WD`, `Mean_P`, `SD_P`, `Skew_P`, `Mean_Tmax`, `SD_Tmax`, `Mean_Tmin`, `SD_Tmin`) in exact order, and the 3 daily outputs (`Precip_Rate`, `Tmax`, `Tmin`).
5. **Run Simulation**:
   Run Monte Carlo realizations. The DLL automatically reseeds reproducibly per realization.

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
