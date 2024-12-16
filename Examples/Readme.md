# Examples
The examples in this folder are to showcase the most fundamental functions of applying stochastic geometry in wireless network, so the selected scenarios and parameters the simplified. For more practical scenarion application, see the folder project

## Interference 
Setup:
* BS density = 10/km^2
* Rayleigh fading with transmission power 1 Watt
* Distance r = 100m

Resutls:
* Average recevived interference power = 
* Computation time = 

## Coverage Probability
Setup:
* BS density = 10/km^2
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Resutls:
* Coverage probability = 
* Computation time = 

Compilation
g++-14 -std=c++20 -I../Integration CoverageProbabilityExample.cpp

## Coverage probability with signals reflected from RISs
Setup:
* BS density = 10/km^2
* RIS cluster = ring of [10, 25]m, 5 RIS per cluster, 
* RIS elements = 200, and Gaussian approximated fading (Rician). 
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Resutls:
* Coverage probability = 
* Computation time = 
* Computation time (Quad core) = 


## Ergodic Rate 
Setup:
* BS density = 10/km^2
* Rayleigh fading with transmission powr 1 Watt
* Distance r = 100m

Resutls:
* Ergodic rate = 
* Computation time = 
* Computation time (Quad core) =  

g++-14 -std=c++20 -I../Integration CoverageProbabilityExample.cpp


Computation time is based on the computer: MacBook Pro 2 GHz Quad core Intel Core i5