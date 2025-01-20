#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include "../LaplaceTransform/PoissonPP.hpp"
#include "../LaplaceTransform/Exponential.hpp"
#include "../LaplaceTransform/ChiSquare.hpp"
#include "../LaplaceTransform/LaplaceTransform.hpp"
#include "../StochasticGeometry/CoverageProbability.hpp"
#include "../Utils/pathloss.hpp"
#include "../Integration/GenzMalik/Cube.hpp"
#include "../Integration/GenzMalik/GM2D.hpp"

int main(){

    // Scenario Setup:
    double noise_power = 1e-13; // Watt 
    double BS_density = 1e-5; //per m^2
    double distance_typical_BS = 80.0; // meter
    double coverage_threshold = 1.0;
    
    double RIS_elements = 200.0;
    double radius_min=10.0;
    double radius_max=25.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);
    double RIS_density = {5.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};

    // Count the timer !
    auto start{std::chrono::steady_clock::now()};

    // Computing the coverage probability in three steps:
    // 1, Define the laplace transform of interference from a point
    auto lt_interference_at_r = [&](std::complex<double> s, double r){
        double interference_power = PathLoss::direct_nlos(r); 
        std::complex<double> return_value = Exponential::eval(s, interference_power);
        return return_value;
    };

    auto lt_reflected_signal_at_y_theta = [&](std::complex<double> s, const double y, const double theta){
        double chi_variance_parameter = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance;
        double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
        double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(distance_typical_BS, y, theta);
        chi_variance_parameter = chi_variance_parameter * pathloss_coefficient;
        std::complex<double> lt_chi = ChiSquare::eval_non_central( - s, non_center_parameter, chi_variance_parameter);
        return lt_chi;
    };
    // 2, Define the laplace transform of interference from a point process
    PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, BS_density, distance_typical_BS);
    PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);

    // 3, Obtaining the coverage probability 
    double coverage_probability = CoverageProbability::calculate_rayleigh_positive(distance_typical_BS, coverage_threshold, \
            noise_power, laplace_transform_interference, laplace_transform_reflected_signals);

    auto end{std::chrono::steady_clock::now()};    
    const std::chrono::duration<double> elapsed_seconds{end - start};

    std::cout << "The coverage probability is " << coverage_probability << std::endl;
    std::cout << "The program runs: " << elapsed_seconds.count() << " seconds." << std::endl;

    return 0;
}