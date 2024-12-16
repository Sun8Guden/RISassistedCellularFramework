#include <cmath>
#include <iostream>
#include "DerivativeCoverageHoleBS.hpp"
#include "DerivativeCoverageHoleRIS.hpp"

#include "../../LaplaceTransform/LaplaceTransform.hpp"
#include "../../Integration/GenzMalik/Cube.hpp"
#include "../../Integration/GenzMalik/GM2D.hpp"
#include "../../LaplaceTransform/PoissonPP.hpp"
#include "../../LaplaceTransform/Exponential.hpp"
#include "../../LaplaceTransform/ChiSquare.hpp"
#include "../../StochasticGeometry/ErgodicRate.hpp"

int main(){

    double noise_power = 1e-13; // Watt 
    // double coverage_threshold = 1.0;
    double RIS_elements = 300.0;
    double radius_min=20.0;
    double radius_max=30.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double reference_BS_density {1.0e-5}; //  10 per km^2
    double reference_coverage_hole {100.0}; // 100 meters

    double direct_link_penalty = 5.0;

    // Basic Laplace transform
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

    double test_BS_density {3.0e-6};
    double RIS_density = {0.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};
    double test_coverage_hole = reference_coverage_hole / std::sqrt(test_BS_density/reference_BS_density);
    double normalization = reference_coverage_hole * std::sqrt(reference_BS_density);

    PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, test_BS_density, test_coverage_hole);
    PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);

    double ergo_rate = ErgodicRate::calculate_rayleigh_positive(test_coverage_hole, laplace_transform_interference, noise_power*direct_link_penalty, laplace_transform_reflected_signals);
    // double ergo_rate = ErgodicRate::calculate_rayleigh(test_coverage_hole, laplace_transform_interference, noise_power*direct_link_penalty);

    std::cout << "The ergodic rate is " << ergo_rate << std::endl;
    return 0;
}