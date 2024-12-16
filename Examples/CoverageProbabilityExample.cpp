#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include "../LaplaceTransform/PoissonPP.hpp"
#include "../LaplaceTransform/Exponential.hpp"
#include "../LaplaceTransform/LaplaceTransform.hpp"
#include "../StochasticGeometry/CoverageProbability.hpp"
#include "../Utils/pathloss.hpp"

int main(){

    // Scenario Setup:
    double noise_power = 1e-13; // Watt 
    double BS_density = 1e-5; //per m^2
    double distance_typical_BS = 100.0; // meter
    double coverage_threshold = 1.0;

    // Count the timer !
    auto start{std::chrono::steady_clock::now()};

    // Computing the coverage probability in three steps:
    // 1, Define the laplace transform of interference from a point
    auto laplace_transform_at_r = [&](std::complex<double> s, double r){
        double interference_power = PathLoss::direct_nlos(r); 
        std::complex<double> return_value = Exponential::eval(s, interference_power);
        return return_value;
    };
    // 2, Define the laplace transform of interference from a point process
    PoissonPP2D_Reduced laplace_transform_interference(laplace_transform_at_r, BS_density, distance_typical_BS);
    // 3, Obtaining the coverage probability 
    double coverage_probability = CoverageProbability::calculate_rayleigh(distance_typical_BS, coverage_threshold, laplace_transform_interference, noise_power);

    auto end{std::chrono::steady_clock::now()};    
    const std::chrono::duration<double> elapsed_seconds{end - start};

    std::cout << "The coverage probability is " << coverage_probability << std::endl;
    std::cout << "The program runs: " << elapsed_seconds.count() << " seconds." << std::endl;

    return 0;
}