# Examples
The examples in this folder are to showcase the most fundamental functions of applying stochastic geometry in wireless network, so the selected scenarios and parameters the simplified. 
Both coverage probability and ergodic rate assessments incorporate scenarios with and without the presence of a Reconfigurable Intelligent Surface (RIS). Specifically, they consider the impact of interference and noise in both cases, while additionally accounting for the influence of reflected signals when the RIS is deployed.
Users can readily adapt the provided examples to investigate more complex scenarios, such as different point processes for BSs and RISs, fading distributions, and customized scenarios. Furthermore, RISs can be replaced by other signal technologies capable of generating auxiliary signals.
For more practical scenarion application, see the folder project.


## Coverage Probability
Setup:
* BS density = 10/km^2
* Transmission powr = 1 Watt
* Rayleigh fading 
* Distance r = 100m

Compilation command (please include the integration directory)
g++-14 -std=c++20 -I../Integration CoverageProbabilityExample.cpp

Resutls:
* Coverage probability = 0.660685
* Computation time = 6e-05 seconds

## Coverage probability with signals reflected from RISs
Setup:
* BS density = 10/km^2
* RIS cluster = ring of [10, 25]m, 5 RIS per cluster, 
* RIS elements = 200, and Gaussian approximated fading (Rician). 
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Compilation: g++-14 -std=c++20 -O3 -I../Integration CoverageProbabilityRISExample.cpp

Resutls:
* Coverage probability = 0.777182
* Computation time = 0.432546 seconds


## Ergodic Rate 
Setup:
* BS density = 10/km^2
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Compilation 
g++-14 -std=c++20 -O3 -I../Integration ErgodicRateRISExample.cpp

Resutls:
* Ergodic rate = 1.11374
* Computation time = 0.00273

## Ergodic rate with signals reflected from RISs

Setup:
* BS density = 10/km^2
* RIS cluster = ring of [10, 25]m, 5 RIS per cluster, 
* RIS elements = 200, and Gaussian approximated fading (Rician). 
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Compilation (please utilize the compiler flag -fopenmp to leverage the parallelism in ergodic rate computations involving reflected signals):
g++-14 -std=c++20 -O3 -fopenmp -I../Integration ErgodicRateRISExample.cpp

* Ergodic rate = 1.27658
* Computation time (Quad core) = 4.22977 seconds

Computation time is based on the computer: MacBook Pro 2 GHz Quad core Intel Core i5