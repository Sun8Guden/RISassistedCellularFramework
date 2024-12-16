#include <cmath>
#include <iostream>
#include "../../LaplaceTransform/LaplaceTransform.hpp"

#include "../../Integration/GenzMalik/Cube.hpp"
#include "../../Integration/GenzMalik/GM2D.hpp"
#include "../../LaplaceTransform/PoissonPP.hpp"
#include "../../LaplaceTransform/Exponential.hpp"
#include "../../LaplaceTransform/ChiSquare.hpp"
#include "../../StochasticGeometry/CoverageProbability.hpp"


int main(){

    double noise_power = 1.0e-13; // Watt 
    double coverage_threshold = 1.0;
    double RIS_elements = 200.0;
    double radius_min=10.0;
    double radius_max=25.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double reference_BS_density {1.004e-5}; //  10 per km^2
    double reference_coverage_hole {80.0}; // 100 meters
    reference_coverage_hole = reference_coverage_hole/std::sqrt(reference_BS_density/1.0e-5);
    double RIS_density = {2.0 / (M_PI *( radius_max*radius_max-radius_min*radius_min) )};

    double direct_link_penalty = 1.0;

    auto lt_interference_at_r = [&](std::complex<double> s, double r){
        double interference_power = PathLoss::direct_nlos(r) * direct_link_penalty; 
        std::complex<double> return_value = Exponential::eval(s, interference_power);
        return return_value;
    };


    auto lt_reflected_signal_at_y_theta = [&](std::complex<double> s, const double y, const double theta){
        double chi_variance_parameter = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
        double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
        double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(reference_coverage_hole, y, theta);
        chi_variance_parameter = chi_variance_parameter * pathloss_coefficient;
        std::complex<double> lt_chi = ChiSquare::eval_non_central( -s, non_center_parameter, chi_variance_parameter);
        return lt_chi;
    };

    PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, reference_BS_density, reference_coverage_hole);
    PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);

    double coverage_probability = CoverageProbability::calculate_rayleigh_positive(reference_coverage_hole, coverage_threshold, noise_power * direct_link_penalty, laplace_transform_interference, laplace_transform_reflected_signals);
    // double coverage_probability = CoverageProbability::calculate_rayleigh(reference_coverage_hole, coverage_threshold, laplace_transform_interference, noise_power);
    std::cout << "The coverage probability under the penalty: " << coverage_probability << std::endl; 

}