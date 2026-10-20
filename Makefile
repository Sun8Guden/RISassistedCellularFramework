# Convenience build for systems without CMake (CMakeLists.txt is the primary build).
#
#   make                       examples, tests, benchmarks and validation programs (serial)
#   make test                  run the regression tests
#   make validate              Monte Carlo validation and the option-pricing cross-check
#   make figures               regenerate docs/data and docs/figures (needs python3 with matplotlib)
#   make paper                 build the research drivers in paper/
#   make OPENMP=1 ...          enable OpenMP (GCC, or Apple Clang with `brew install libomp`)
#   make clean
#
# Binaries go to build/, mirroring the source tree.

CXXFLAGS ?= -std=c++20 -O3 -Wall
CPPFLAGS += -Iinclude -Ithird_party

ifeq ($(OPENMP),1)
  ifneq (,$(findstring Apple,$(shell $(CXX) --version 2>/dev/null)))
    LIBOMP ?= $(shell brew --prefix libomp 2>/dev/null)
    CXXFLAGS += -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    LDLIBS += -L$(LIBOMP)/lib -lomp
  else
    CXXFLAGS += -fopenmp
  endif
endif

BUILD := build
DATA := docs/data

binaries = $(patsubst %.cpp,$(BUILD)/%,$(wildcard $(1)))
EXAMPLES := $(call binaries,examples/*.cpp)
TESTS := $(call binaries,tests/*.cpp)
BENCHMARKS := $(call binaries,benchmarks/*.cpp)
VALIDATION := $(call binaries,validation/*.cpp)
PAPER := $(call binaries,paper/*/*.cpp paper/*/simulation_*/*.cpp)
HEADERS := $(wildcard include/risgeo/*.hpp include/risgeo/*/*.hpp include/risgeo/*/*/*.hpp examples/*/*.hpp \
                      paper/*/*.hpp paper/*/*/*.hpp)

.PHONY: all test validate figures paper clean

all: $(EXAMPLES) $(TESTS) $(BENCHMARKS) $(VALIDATION)

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done

validate: $(BUILD)/validation/monte_carlo $(BUILD)/examples/OptionPricingExample
	./$(BUILD)/validation/monte_carlo
	./$(BUILD)/examples/OptionPricingExample

figures: $(BUILD)/validation/monte_carlo $(BUILD)/benchmarks/cubature_benchmark $(BUILD)/examples/OptionPricingExample
	@mkdir -p $(DATA)
	./$(BUILD)/validation/monte_carlo 4000000 $(DATA)
	./$(BUILD)/benchmarks/cubature_benchmark $(DATA)/cubature_convergence.csv
	./$(BUILD)/examples/OptionPricingExample $(DATA)
	python3 scripts/make_figures.py $(DATA) docs/figures

paper: $(PAPER)

$(BUILD)/%: %.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -rf $(BUILD)
