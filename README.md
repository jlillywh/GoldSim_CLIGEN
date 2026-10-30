# GoldSim - CLIGEN Integration

[![Windows 64-bit](https://img.shields.io/badge/platform-Windows%2064--bit-blue.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`CligenWeather.dll` is a compiled C++ implementation of the **CLIGEN** stochastic weather generator that interfaces directly with [GoldSim](https://www.goldsim.com/Web/Home/)'s External (DLL) element to simulate daily weather data (precipitation, $T_{\max}$, and $T_{\min}$) during Monte Carlo simulations.

---

## About CLIGEN

**CLIGEN** (Climate Generator) is a stochastic weather generator developed by the **USDA Agricultural Research Service (ARS)** originally for the Water Erosion Prediction Project (WEPP) and widely used in hydrological, agricultural, and environmental modeling (Nicks et al., 1995; Flanagan & Nearing, 1995). Rather than replaying a single historical weather record, CLIGEN uses monthly statistical parameters derived from long-term meteorological station records to generate synthetic daily weather sequences of arbitrary length across multiple Monte Carlo realizations.

Key statistical components of the generator include:

1. **Precipitation Occurrence (First-Order Two-State Markov Chain)**:
   Daily wet and dry state transitions are modeled using two monthly conditional probabilities: $P(W|W)$, the probability of a wet day following a wet day, and $P(W|D)$, the probability of a wet day following a dry day. This structure preserves the seasonal frequency of precipitation as well as the distribution of wet-spell and dry-spell durations.

2. **Precipitation Depth (Skewed Pearson Type III Distribution)**:
   On simulated wet days, daily precipitation depth is sampled from a three-parameter **Pearson Type III** (skewed Gamma-family) distribution using the Wilson–Hilferty transformation, parameterized by each month's wet-day mean (`Mean_P`), standard deviation (`SD_P`), and skewness coefficient (`Skew_P`). Incorporating historical skewness reproduces both frequent light rainfall events and extreme right-tail storm depths.

3. **Daily Air Temperature ($T_{\max}$ and $T_{\min}$)**:
   Daily maximum and minimum temperatures are generated from monthly means (`Mean_Tmax`, `Mean_Tmin`) and standard deviations (`SD_Tmax`, `SD_Tmin`) using correlated normal variates with lag-1 day-to-day persistence ($\rho = 0.60$) and same-day cross-correlation ($r = 0.65$). Maximum temperature is conditioned on the day's wet/dry state to reflect reduced solar heating and evaporative cooling on rainy days while preserving the unconditional monthly mean.

4. **Continuous Seasonal Interpolation**:
   To avoid step changes at month boundaries, `CligenWeather.dll` fits a periodic cubic spline through the 12 monthly values of each input parameter vector and evaluates the smooth seasonal curve at the current day of the year on each simulation timestep.

---

## GoldSim Integration Steps

### 1. Drop DLL into Model Directory
Place `CligenWeather.dll` in the directory containing your GoldSim (`.gsm`) model. No installation or external runtime environment is required.
* **64-bit GoldSim**: Use `CligenWeather.dll` (or `CligenWeather_x64.dll`) in-process.
* **32-bit GoldSim**: Check **"Run in separate process"** on the External element's Definition tab to avoid Windows Error 193 (`ERROR_BAD_EXE_FORMAT`), or use the 32-bit build `CligenWeather_win32.dll` (renamed to `CligenWeather.dll`).

### 2. Paste Parameters into Data Elements
Define 12-item monthly vector parameters (Data elements indexed by `Months`, Jan–Dec) in GoldSim to supply the 9 monthly parameter vectors (`P_WW`, `P_WD`, `Mean_P`, `SD_P`, `Skew_P`, `Mean_Tmax`, `SD_Tmax`, `Mean_Tmin`, `SD_Tmin`) to the External element.

### 3. Link External Element & Run
In GoldSim, insert an External element (**Insert Element → Functions → External (DLL)**):
* **File Name**: `CligenWeather.dll`
* **Function Name**: `CligenWeather` (or `cligen_weather_step` / `cligen_weather_arrays`)
* **Calling Convention**: `__cdecl`

Link the 4 run-property scalars, 9 monthly 12-item vectors, and 3 daily outputs (`Precip_Rate`, `Tmax`, `Tmin`).

---

## Technical Interface Specification

### DLL Input Arguments (`GoldSim → DLL`)
Values GoldSim sends to the DLL at each simulation step (4 scalars + 9 monthly 12-item vectors):

| Input Argument | Type / Unit | Description |
| :--- | :--- | :--- |
| `Simulation_Time` | `Scalar (day)` | Current elapsed simulation time (link to GoldSim `ETIME` or `DateTime`) |
| `Time_Step` | `Scalar (day)` | Model time step (link to GoldSim `DT`, e.g. 1.0 day) |
| `Realization` | `Scalar (-)` | Current Monte Carlo realization (link to GoldSim `Realization`) |
| `Units_Flag` | `Scalar (0 or 1)` | `0` = Metric (mm, &deg;C), `1` = US Customary (in, &deg;F) |
| `P_WW` | `Vector[12] (-)` | Monthly wet-day probability given yesterday was wet (12 items, Jan–Dec) |
| `P_WD` | `Vector[12] (-)` | Monthly wet-day probability given yesterday was dry (12 items, Jan–Dec) |
| `Mean_P` | `Vector[12] (in or mm)` | Monthly mean daily precipitation depth on wet days (12 items, Jan–Dec) |
| `SD_P` | `Vector[12] (in or mm)` | Monthly daily precipitation standard deviation (12 items, Jan–Dec) |
| `Skew_P` | `Vector[12] (-)` | Monthly daily precipitation skewness coefficient (12 items, Jan–Dec) |
| `Mean_Tmax` | `Vector[12] (°F or °C)` | Monthly mean daily maximum air temperature (12 items, Jan–Dec) |
| `SD_Tmax` | `Vector[12] (Fdeg or Cdeg)` | Monthly daily maximum temperature standard deviation (12 items, Jan–Dec) |
| `Mean_Tmin` | `Vector[12] (°F or °C)` | Monthly mean daily minimum air temperature (12 items, Jan–Dec) |
| `SD_Tmin` | `Vector[12] (Fdeg or Cdeg)` | Monthly daily minimum temperature standard deviation (12 items, Jan–Dec) |

### DLL Outputs Returned to GoldSim (`DLL → GoldSim`)
Values updated by the DLL on each daily simulation timestep (3 values):

| Output Parameter | Units | Description |
| :--- | :--- | :--- |
| `Precip_Rate` | `mm/day or in/day` | Simulated stochastic precipitation rate for the current timestep (`prcp / DT`) |
| `Tmax` | `°C or °F` | Simulated daily maximum air temperature |
| `Tmin` | `°C or °F` | Simulated daily minimum air temperature |

---

## Logging & Diagnostics

By default, `CligenWeather.dll` operates in **`ERROR`** logging mode with zero disk I/O during normal timesteps.
* **Default (`ERROR`)**: Only critical errors (such as null pointers) are written to `cligen_debug.log`.
* **Enable Debug Logging (`DEBUG`)**: Set the Windows environment variable `CLIGEN_LOG_LEVEL=DEBUG` or create an empty file named `cligen_debug.txt` in the same directory as your `.gsm` model to write step-by-step diagnostics to `cligen_debug.log`.
* **Disable All Logging (`OFF`)**: Set `CLIGEN_LOG_LEVEL=OFF`.

---

## Building from Source

Pre-compiled binaries (`CligenWeather.dll`, `CligenWeather_x64.dll`, and `CligenWeather_win32.dll`) are included in the repository. To compile from source:

### Windows (MSVC or MinGW)
```cmd
build.bat
```
Or manually:
* **MinGW-w64 (`g++`)**:
  ```cmd
  g++ -O3 -std=c++17 -Wall -shared goldsim_external.cpp -o CligenWeather.dll -static -Wl,--add-stdcall-alias
  ```
* **Microsoft Visual C++ (`cl`)**:
  ```cmd
  cl /O2 /std:c++17 /EHsc /LD goldsim_external.cpp /Fe:CligenWeather.dll
  ```

### Linux (Cross-Compilation via MinGW-w64)
```sh
make dll
```

---

## References & Further Reading

* Nicks, A. D., Lane, L. J., & Gander, G. A. (1995). *Weather Generator*. Chapter 2 in **USDA–Water Erosion Prediction Project: Hillslope Profile and Watershed Model Documentation**, NSERL Report No. 10, USDA-ARS National Soil Erosion Research Laboratory, West Lafayette, IN.
* Flanagan, D. C., & Nearing, M. A. (Eds.) (1995). *USDA–Water Erosion Prediction Project: Hillslope Profile and Watershed Model Documentation*. NSERL Report No. 10, USDA-ARS National Soil Erosion Research Laboratory.
* Zhang, X. C., & Garbrecht, J. D. (2003). *Evaluation of CLIGEN precipitation parameters and their implication on WEPP runoff and erosion prediction*. Transactions of the ASAE, 46(2), 311–320.
* Kou, X., Ge, J., Wang, Y., & Zhang, C. (2007). *Validation of the weather generator CLIGEN with daily precipitation data from the Loess Plateau, China*. Journal of Hydrology, 347(3–4), 347–357.
* Meyer, C. R., Renschler, C. S., & Vining, R. C. (2008). *Implementing quality control on a random number stream to improve a stochastic weather generator*. Hydrological Processes, 22(8), 1069–1079.
* USDA-ARS National Soil Erosion Research Laboratory: [CLIGEN Documentation & Overview](https://www.ars.usda.gov/midwest-area/west-lafayette-in/national-soil-erosion-research/docs/wepp/cligen/).

---

## License

MIT License &mdash; Free for academic, personal, and commercial engineering use.
Developed by [Lillywhite Consulting](https://lillywhiteconsulting.com).
