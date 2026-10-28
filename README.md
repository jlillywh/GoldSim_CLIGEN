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
* **Exported Function**: `cligen_weather_step`
* **Calling Convention**: `__cdecl` (Standard C)

### Inputs (`inargs` &mdash; 13 Scalars)
| Index | Argument Name | GoldSim Type | Typical Unit | Description |
| :---: | :--- | :--- | :--- | :--- |
| `[0]` | `Simulation_Time` | Scalar Value | `day` | Elapsed simulation time (GoldSim `ETIME`) |
| `[1]` | `Time_Step` | Scalar Value | `day` | Model timestep size (GoldSim `DT`, e.g. 1.0 day) |
| `[2]` | `Realization` | Scalar Value | `-` | Current Monte Carlo realization (GoldSim `Realization`) |
| `[3]` | `Units_Flag` | Scalar Value | `-` | `0` = Metric (mm, &deg;C), `1` = US Customary (in, &deg;F) |
| `[4]` | `P_WW` | Scalar Value | `-` | Wet-day probability given yesterday was wet $P(W\|W)$ |
| `[5]` | `P_WD` | Scalar Value | `-` | Wet-day probability given yesterday was dry $P(W\|D)$ |
| `[6]` | `Mean_P` | Scalar Value | `in` or `mm` | Daily mean precipitation depth on wet days |
| `[7]` | `SD_P` | Scalar Value | `in` or `mm` | Daily precipitation standard deviation |
| `[8]` | `Skew_P` | Scalar Value | `-` | Daily precipitation skewness coefficient |
| `[9]` | `Mean_Tmax` | Scalar Value | `&deg;F` or `&deg;C` | Daily mean maximum air temperature |
| `[10]` | `SD_Tmax` | Scalar Value | `Fdeg` or `Cdeg` | Daily maximum temperature standard deviation |
| `[11]` | `Mean_Tmin` | Scalar Value | `&deg;F` or `&deg;C` | Daily mean minimum air temperature |
| `[12]` | `SD_Tmin` | Scalar Value | `Fdeg` or `Cdeg` | Daily minimum temperature standard deviation |

### Outputs (`outargs` &mdash; 3 Scalars)
| Index | Output Parameter | Units | Description |
| :---: | :--- | :--- | :--- |
| `[0]` | `Precip_Rate` | `in/day` or `mm/day` | Stochastic precipitation rate for the current timestep (`prcp / DT`) |
| `[1]` | `Tmax` | `&deg;F` or `&deg;C` | Daily maximum air temperature |
| `[2]` | `Tmin` | `&deg;F` or `&deg;C` | Daily minimum air temperature |

---

## GoldSim Setup Guide

1. **Place DLL**: Copy `CligenWeather.dll` into the directory containing your GoldSim model file (`.gsm`).
2. **Define Scalar Parameters**:
   Create scalar elements in GoldSim (such as Data elements, Expression elements, or Lookup Tables) for each of the 9 climate parameters:
   * `P_WW` &mdash; Daily wet-given-wet probability (scalar fraction)
   * `P_WD` &mdash; Daily wet-given-dry probability (scalar fraction)
   * `Mean_P` &mdash; Daily wet-day mean precipitation depth (scalar length: `in` or `mm`)
   * `SD_P` &mdash; Daily precipitation standard deviation (scalar length: `in` or `mm`)
   * `Skew_P` &mdash; Daily precipitation skewness coefficient (scalar fraction)
   * `Mean_Tmax` &mdash; Daily mean maximum air temperature (scalar temperature: `&deg;F` or `&deg;C`)
   * `SD_Tmax` &mdash; Daily maximum temperature standard deviation (scalar temp difference: `Fdeg` or `Cdeg`)
   * `Mean_Tmin` &mdash; Daily mean minimum air temperature (scalar temperature: `&deg;F` or `&deg;C`)
   * `SD_Tmin` &mdash; Daily minimum temperature standard deviation (scalar temp difference: `Fdeg` or `Cdeg`)

   *(If your model uses monthly values, you can look them up by calendar month or day so the current scalar value for today is fed to the element.)*

3. **Insert External Element**:
   * Right-click canvas &rarr; **Insert Element &rarr; Functions &rarr; External (DLL)**.
   * File Name: `CligenWeather.dll`
   * Function: `cligen_weather_step`
   * Calling Convention: `__cdecl`
4. **Link Arguments**:
   Link the 13 scalar inputs in order. Every single argument in the External element interface is a pure scalar value.
5. **Run**:
   Run Monte Carlo simulations across any number of realizations. The DLL automatically reseeds reproducibly per realization.

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
