# Makefile — CLIGEN pure-C++ GoldSim External element
# ---------------------------------------------------------------------------
# Targets:
#   make test     build the standalone verification harness (no Python)
#   make dll      build the GoldSim External DLL (Windows / MinGW cross-build)
#   make so       build a Linux shared object (local testing)
#   make run      build + run the 30-year statistical sanity check
#   make clean    remove build artifacts
#
# The DLL/so export the single GoldSim entry point CligenWeather().
# The test harness links the engine statically and drives it exactly as
# GoldSim would, day by day.
# ---------------------------------------------------------------------------

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall
ENGINE    = goldsim_external.cpp
HEADER    = cligen_engine.hpp

.PHONY: test run dll so clean

test: cligen_test

cligen_test: $(ENGINE) test_harness.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) $(ENGINE) test_harness.cpp -o $@

# Linux shared object (quick local load test).
so: cligen_engine.so
cligen_engine.so: $(ENGINE) $(HEADER)
	$(CXX) $(CXXFLAGS) -shared -fPIC $(ENGINE) -o $@

# Windows DLLs for GoldSim:
#   make dll32   build 32-bit DLL for standard 32-bit GoldSim (in-process, fastest)
#   make dll64   build 64-bit DLL for 64-bit GoldSim
#   make dll     build both
CXX32 ?= i686-w64-mingw32-g++
CXX64 ?= x86_64-w64-mingw32-g++

dll: dll32 dll64

dll32: CligenWeather_win32.dll
CligenWeather_win32.dll: $(ENGINE) $(HEADER)
	$(CXX32) -O3 -std=c++17 -Wall -shared $(ENGINE) -o $@ -static -Wl,--add-stdcall-alias
	cp $@ CligenWeather.dll

dll64: CligenWeather_x64.dll
CligenWeather_x64.dll: $(ENGINE) $(HEADER)
	$(CXX64) -O3 -std=c++17 -Wall -shared $(ENGINE) -o $@ -static -Wl,--add-stdcall-alias

run: cligen_test
	mkdir -p ../../../bridge_output
	./cligen_test --years 30 --start-year 1995 --seed 1 \
	    --out ../../../bridge_output/cpp_sim_output.csv

clean:
	rm -f cligen_test cligen_engine.so CligenWeather.dll CligenWeather_win32.dll CligenWeather_x64.dll
