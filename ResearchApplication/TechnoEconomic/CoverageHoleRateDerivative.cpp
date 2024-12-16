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

int main(){

    double noise_power = 1e-13; // Watt 
    double coverage_threshold = 1.0;
    double RIS_elements = 300.0;
    double radius_min=10.0;
    double radius_max=25.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double reference_BS_density {1.0e-5}; //  10 per km^2
    double reference_coverage_hole {80.0}; // 100 meters
    double RIS_density = {0.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};

    double direct_link_penalty = 1.0;

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

    // The derivative of the reflected signal part! 




    // 2, Define the laplace transform of interference from a point process

    double test_BS_density {1.0e-5};
    double test_coverage_hole = reference_coverage_hole / std::sqrt(test_BS_density/reference_BS_density);
    double normalization = reference_coverage_hole * std::sqrt(reference_BS_density);



    auto lt_reflected_signal_derivative_s_r = [&](std::complex<double> s, const double y, const double theta){
        double chi_variance_parameter_normalize = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
        double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
        double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(test_coverage_hole, y, theta);
        
        double chi_variance_parameter = chi_variance_parameter_normalize * pathloss_coefficient;
        std::complex<double> denominator =  1.0 - 2.0 * chi_variance_parameter * s;
        std::complex<double> lt_chi_D1 = std::exp( s * non_center_parameter * chi_variance_parameter / denominator )/( std::sqrt(denominator) * denominator * denominator);
        lt_chi_D1 = lt_chi_D1 * (1.0 + non_center_parameter - 2.0 * s * chi_variance_parameter)*\
            (chi_variance_parameter * Derivative::partial_argument_BS(test_coverage_hole, test_BS_density, normalization) + \
            s * chi_variance_parameter_normalize * PathLoss::direct_los(y) * PathLoss::reflected_los_derivative(test_coverage_hole, y, theta) * Derivative::partial_CH_center_BS(test_BS_density, normalization));
        
        return lt_chi_D1; 
    };

    auto lt_reflected_signal_derivative_r = [&](std::complex<double> s, const double y, const double theta){
        double chi_variance_parameter_normalize = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
        double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
        double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(test_coverage_hole, y, theta);
        
        double chi_variance_parameter = chi_variance_parameter_normalize * pathloss_coefficient;
        std::complex<double> denominator =  1.0 - 2.0 * chi_variance_parameter * s;
        std::complex<double> lt_chi_D1 = std::exp( s * non_center_parameter * chi_variance_parameter / denominator )/( std::sqrt(denominator) * denominator * denominator);
        lt_chi_D1 = lt_chi_D1 * (1.0 + non_center_parameter - 2.0*s* chi_variance_parameter)*\
            (s * chi_variance_parameter_normalize * PathLoss::direct_los(y) * \
            PathLoss::reflected_los_derivative(test_coverage_hole, y, theta)*Derivative::partial_CH_center_BS(test_BS_density, normalization));
    
        return lt_chi_D1; 
    };

    PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, test_BS_density, test_coverage_hole);
    PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);
    PoissonPP2D_partial laplace_transform_reflected_signals_partial_s_r(lt_reflected_signal_derivative_s_r, RIS_density, cluster_cube);
    PoissonPP2D_partial laplace_transform_reflected_signals_partial_r(lt_reflected_signal_derivative_r, RIS_density, cluster_cube);

    // For any given BS density, the distance to the associated BS is normalized to guarantee that the relative distance between BS and RIS remains the constant.

    double derivative_BS = Derivative::get_BS_derivative_rate(test_coverage_hole, test_BS_density, normalization, \
            noise_power * direct_link_penalty, laplace_transform_interference, laplace_transform_reflected_signals, \
            laplace_transform_reflected_signals_partial_s_r, laplace_transform_reflected_signals_partial_r);
    // double derivative_BS_no_RIS = Derivative::get_BS_coverage_no_RIS(test_coverage_hole, coverage_threshold, noise_power, laplace_transform_interference, test_BS_density, normalization);
    double derivatove_RIS = Derivative::get_RIS_derivative_rate(test_coverage_hole, noise_power * direct_link_penalty, \
            laplace_transform_interference, laplace_transform_reflected_signals); 
    
        // double derivatove_RIS = Derivative::get_RIS_derivative(test_coverage_hole, 1.0, noise_power, \
        //     laplace_transform_interference, laplace_transform_reflected_signals); 
    std::cout << "Ergodic Rate is " << std::endl;
    std::cout << "Derivative w.r.t. BS is " << derivative_BS * 1.0e-6 << std::endl;
    // std::cout << "Derivative w.r.t. BS (direct) is " << derivative_BS_no_RIS * 1.0e-6 << std::endl;

    derivatove_RIS = derivatove_RIS / (M_PI*(radius_max*radius_max-radius_min*radius_min)); // This is normalized for adding one RIS per cluster. 
    std::cout << "Derivative w.r.t. RIS is " << derivatove_RIS << std::endl;

    return 0; 
}