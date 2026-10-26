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

# Windows DLL for GoldSim. Use a MinGW toolchain, e.g.:
#   make dll CXX=x86_64-w64-mingw32-g++
# or build natively on Windows with the same source.
dll: CligenWeather.dll
CligenWeather.dll: $(ENGINE) $(HEADER)
	$(CXX) $(CXXFLAGS) -shared $(ENGINE) -o $@ -static -Wl,--add-stdcall-alias

run: cligen_test
	mkdir -p ../../../bridge_output
	./cligen_test --years 30 --start-year 1995 --seed 1 \
	    --out ../../../bridge_output/cpp_sim_output.csv

clean:
	rm -f cligen_test cligen_engine.so CligenWeather.dll
