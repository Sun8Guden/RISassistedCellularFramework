#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>

#include "../DerivativeCoverageHoleBS.hpp"
#include "../DerivativeCoverageHoleRIS.hpp"

#include "../../../LaplaceTransform/LaplaceTransform.hpp"
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"
#include "../../../LaplaceTransform/PoissonPP.hpp"
#include "../../../LaplaceTransform/Exponential.hpp"
#include "../../../LaplaceTransform/ChiSquare.hpp"


#define _USE_MATH_DEFINES // for using PI in C++

#include "../../../Utils/json.hpp"
using json = nlohmann::json;

// This simulation investigate the impact of direct link penalty on the decision selection.

int main(){ 

    double noise_power = 1e-13; // Watt 
    double coverage_threshold = 1.0;
    double RIS_elements = 600.0;
    double radius_min=20.0;
    double radius_max=30.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double reference_BS_density {1.0e-5}; //  10 per km^2
    double reference_coverage_hole {80.0}; // 100 meters
    double normalization = reference_coverage_hole * std::sqrt(reference_BS_density);

    double RIS_density = {3.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};

    std::vector<double> rate_gain_BS;
    std::vector<double> rate_gain_RIS;


    for (int i=0; i<11; i++){
        double direct_link_penalty = std::pow(10.0, double(i)/10.0); // 
        std::cout << "Round " << i+1 << ": Penalty is " << direct_link_penalty << std::endl;
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

        auto lt_reflected_signal_derivative_s_r = [&](std::complex<double> s, const double y, const double theta){
            double chi_variance_parameter_normalize = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
            double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
            double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(reference_coverage_hole, y, theta);
            
            double chi_variance_parameter = chi_variance_parameter_normalize * pathloss_coefficient;
            std::complex<double> denominator =  1.0 - 2.0 * chi_variance_parameter * s;
            std::complex<double> lt_chi_D1 = std::exp( s * non_center_parameter * chi_variance_parameter / denominator )/( std::sqrt(denominator) * denominator * denominator);
            lt_chi_D1 = lt_chi_D1 * (1.0 + non_center_parameter - 2.0 * s * chi_variance_parameter)*\
                (chi_variance_parameter * Derivative::partial_argument_BS(reference_coverage_hole, reference_BS_density, normalization) + \
                s * chi_variance_parameter_normalize * PathLoss::direct_los(y) * \
                PathLoss::reflected_los_derivative(reference_coverage_hole, y, theta) * Derivative::partial_CH_center_BS(reference_BS_density, normalization));
            
            return lt_chi_D1; 
        };


        auto lt_reflected_signal_derivative_r = [&](std::complex<double> s, const double y, const double theta){
            double chi_variance_parameter_normalize = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
            double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
            double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(reference_coverage_hole, y, theta);
            
            double chi_variance_parameter = chi_variance_parameter_normalize * pathloss_coefficient;
            std::complex<double> denominator =  1.0 - 2.0 * chi_variance_parameter * s;
            std::complex<double> lt_chi_D1 = std::exp( s * non_center_parameter * chi_variance_parameter / denominator )/( std::sqrt(denominator) * denominator * denominator);
            lt_chi_D1 = lt_chi_D1 * (1.0 + non_center_parameter - 2.0*s* chi_variance_parameter)*\
                (s * chi_variance_parameter_normalize * PathLoss::direct_los(y) * \
                PathLoss::reflected_los_derivative(reference_coverage_hole, y, theta)*Derivative::partial_CH_center_BS(reference_BS_density, normalization));
        
            return lt_chi_D1; 
        };

        PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, reference_BS_density, reference_coverage_hole);
        PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);
        PoissonPP2D_partial laplace_transform_reflected_signals_partial_s_r(lt_reflected_signal_derivative_s_r, RIS_density, cluster_cube);
        PoissonPP2D_partial laplace_transform_reflected_signals_partial_r(lt_reflected_signal_derivative_r, RIS_density, cluster_cube);

        // For any given BS density, the distance to the associated BS is normalized to guarantee that the relative distance between BS and RIS remains the constant.

        double derivative_BS = Derivative::get_BS_derivative_rate(reference_coverage_hole, reference_BS_density, normalization, \
                noise_power * direct_link_penalty, laplace_transform_interference, laplace_transform_reflected_signals, \
                laplace_transform_reflected_signals_partial_s_r, laplace_transform_reflected_signals_partial_r);
        // double derivative_BS_no_RIS = Derivative::get_BS_coverage_no_RIS(test_coverage_hole, coverage_threshold, noise_power, laplace_transform_interference, test_BS_density, normalization);
        double derivative_RIS = Derivative::get_RIS_derivative_rate(reference_coverage_hole, noise_power * direct_link_penalty, \
                laplace_transform_interference, laplace_transform_reflected_signals); 

        double normalized_gain_BS = derivative_BS * 1e-6;
        double normailzed_gain_RIS = derivative_RIS / (M_PI*(radius_max*radius_max-radius_min*radius_min)) / (reference_BS_density * 1e6);

        rate_gain_BS.push_back(normalized_gain_BS);
        rate_gain_RIS.push_back(normailzed_gain_RIS);

    }

    json gain_BS_output(rate_gain_BS);
    json gain_RIS_output(rate_gain_RIS);

    std::ofstream out_file( "output.json" );
    out_file << "BS_improvement = " << gain_BS_output << ";" <<std::endl;
    out_file << "RIS_improvement  = " << gain_RIS_output << ";" <<std::endl;



    return 0;
}