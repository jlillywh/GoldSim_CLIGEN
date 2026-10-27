# GoldSim - CLIGEN Integration

[![Windows 64-bit](https://img.shields.io/badge/platform-Windows%2064--bit-blue.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://github.com/jlillywh/GoldSim_CLIGEN)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[GoldSim](https://www.goldsim.com/Web/Home/) DLL `CligenWeather.dll` (compiled C++ Windows 64-bit). Links directly to GoldSim's External Element to run stochastic simulations of daily and hourly weather data.

---

## Key Highlights

* **Zero External Dependencies**: Standard C++17 only. No Python shims, no DLL-to-Python marshaling, no external runtime files.
* **Step-by-Step Scalar Interface**: GoldSim passes 13 scalar inputs at each daily timestep and receives 3 daily scalar outputs (`Precipitation`, `Tmax`, `Tmin`).
* **Clean Unit Handling in GoldSim**: Because inputs are individual scalars, each parameter carries its native GoldSim dimension (`F` for temperature, `Fdeg` for temperature standard deviation, `in` or `mm` for precipitation) with **zero array dimensional conflicts**.
* **Dynamic / Transient Climate Shifts**: Modeler can apply time-varying trends (e.g. warming offsets or shifting precipitation regimes) directly in GoldSim before passing daily parameters to the DLL.
* **CLIGEN Physical Formulation**:
  * Two-state first-order Markov chain for wet/dry day occurrence ($P(W|W)$ and $P(W|D)$).
  * Pearson Type III skew distribution for daily precipitation depth via Wilson-Hilferty transformation.
  * Bivariate normal distribution with cross-correlation ($r = 0.65$) and $\text{AR}(1)$ daily persistence ($\rho = 0.60$) for maximum and minimum temperatures.
  * Evaporative cooling temperature depression on wet days.

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
| `[0]` | `Precipitation` | `in` or `mm` | Daily stochastic precipitation depth |
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
