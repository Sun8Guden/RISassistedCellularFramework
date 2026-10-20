# Examples
The examples in this folder are to showcase the most fundamental functions of applying stochastic geometry in wireless network, so the selected scenarios and parameters are simplified. 
Both coverage probability and ergodic rate assessments incorporate scenarios with and without the presence of a Reconfigurable Intelligent Surface (RIS). Specifically, they consider the impact of interference and noise in both cases, while additionally accounting for the influence of reflected signals when the RIS is deployed.
Users can readily adapt the provided examples to investigate more complex scenarios, such as different point processes for BSs and RISs, fading distributions, and customized scenarios. Furthermore, RISs can be replaced by other signal technologies capable of generating auxiliary signals.
For more practical scenario applications, see the folder `paper/`.

Build from the repository root with `cmake -S . -B build && cmake --build build` (or `make`); binaries land in `build/examples/`. OpenMP parallelizes the ergodic rate with reflected signals when it is available. The analytical results are validated against Monte Carlo simulation in `validation/`.

| Program | What it shows |
| --- | --- |
| `QuickStart.cpp` | the shortest complete use of the library (quoted in the top-level README) |
| `CoverageProbabilityExample.cpp`, `CoverageProbabilityRISExample.cpp` | coverage probability without and with RISs |
| `ErgodicRateExample.cpp`, `ErgodicRateRISExample.cpp` | ergodic rate without and with RISs |
| `OptionPricingExample.cpp` | the same Gil-Pelaez inversion pricing Black-Scholes and Heston options (`finance/pricing.hpp`) |

Computation times below were measured on a MacBook Pro 2 GHz Quad core Intel Core i5 (see the top-level README for Apple M4 timings).

## Coverage Probability
Setup:
* BS density = 10/km^2
* Transmission power = 1 Watt
* Rayleigh fading 
* Distance r = 100m

Results:
* Coverage probability = 0.660685
* Computation time = 6e-05 seconds

## Coverage probability with signals reflected from RISs
Setup:
* BS density = 10/km^2
* RIS cluster = ring of [10, 25]m, 5 RIS per cluster, 
* RIS elements = 200, and Gaussian approximated fading (Rician). 
* Rayleigh fading with transmission power 1 Watt
* Distance r = 100m

Results:
* Coverage probability = 0.777182
* Computation time = 0.432546 seconds


## Ergodic Rate 
Setup:
* BS density = 10/km^2
* Rayleigh fading with transmission power 1 Watt
* Distance r = 100m

Results:
* Ergodic rate = 1.11374 nat/s/Hz
* Computation time = 0.00273 seconds

## Ergodic rate with signals reflected from RISs

Setup:
* BS density = 10/km^2
* RIS cluster = ring of [10, 25]m, 5 RIS per cluster, 
* RIS elements = 200, and Gaussian approximated fading (Rician). 
* Rayleigh fading with transmission power 1 Watt
* Distance r = 100m

Results:
* Ergodic rate = 1.27658 nat/s/Hz
* Computation time (Quad core, OpenMP) = 4.22977 seconds
