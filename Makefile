# ===========================================================================
# CLIGEN C++ Engine & GoldSim External DLL Build Configuration
# ===========================================================================

CXX64 = x86_64-w64-mingw32-g++
HOST_CXX = g++

CXXFLAGS = -O3 -std=c++17 -Wall

# Core source files
ENGINE = goldsim_external.cpp
TEST_SRC = test_harness.cpp goldsim_external.cpp

.PHONY: all dll test clean

all: dll test

# Build the single 64-bit Windows DLL for GoldSim External Element (Time Series output)
dll:
	@echo "Building 64-bit Windows DLL (CligenWeather.dll)..."
	$(CXX64) $(CXXFLAGS) -shared $(ENGINE) -o CligenWeather.dll -static -Wl,--add-stdcall-alias
	@echo "Done: CligenWeather.dll"

# Build and run the native Linux test harness to verify DLL entry points
test:
	@echo "Building and running native Linux test harness..."
	$(HOST_CXX) $(CXXFLAGS) $(TEST_SRC) -o test_cligen
	./test_cligen

clean:
	rm -f CligenWeather.dll test_cligen cligen_debug.log
