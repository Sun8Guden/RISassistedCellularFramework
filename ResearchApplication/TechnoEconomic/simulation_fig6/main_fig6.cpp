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
#include "../../../StochasticGeometry/ErgodicRate.hpp"



#define _USE_MATH_DEFINES // for using PI in C++

#include "../../../Utils/json.hpp"
using json = nlohmann::json;

int main(){
    double cost_ratio = 10.0;
    double BS_increment = 1.0e-6;

    double noise_power = 1e-13; // Watt 
    // double coverage_threshold = 1.0;
    double RIS_elements = 600.0;
    double radius_min=20.0;
    double radius_max=30.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double direct_link_penalty = 1.0;
    double reference_BS_density {1.0e-5}; //  10 per km^2
    double reference_coverage_hole {80.0}; // 100 meters
    double normalization = reference_coverage_hole * std::sqrt(reference_BS_density);


    // double RIS_density = {0.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};

    auto lt_interference_at_r = [&](std::complex<double> s, double r){
        double interference_power = PathLoss::direct_nlos(r) * direct_link_penalty; 
        std::complex<double> return_value = Exponential::eval(s, interference_power);
        // std::cout<< "Value " << return_value << std::endl;
        return return_value;
    };



    std::vector<double> ergo_rate_vec;
    std::vector<double> BS_gain_vec;
    std::vector<double> RIS_gain_vec;
    std::vector<bool>   decision_vec;
    std::vector<double> BS_density_vec;
    std::vector<double> RIS_density_vec;

    double BS_density = 5.0e-6;
    double RIS_density = {0.0 / (M_PI*(radius_max*radius_max-radius_min*radius_min))};


    for (int investment_round = 0; investment_round < 10; investment_round++){

        std::cout << "Investment Round: " << investment_round << " Round!\n";
        std::cout << "BS density is " << BS_density << ", RIS density is " << RIS_density << std::endl;
        BS_density_vec.push_back(BS_density);
        RIS_density_vec.push_back(RIS_density);

        double test_coverage_hole = reference_coverage_hole / std::sqrt(BS_density/reference_BS_density);

        // std::cout << "test_coverage_hole "<< test_coverage_hole << std::endl;

        auto lt_reflected_signal_at_y_theta = [&](std::complex<double> s, const double y, const double theta){
                double chi_variance_parameter = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
                double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
                double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(test_coverage_hole, y, theta);
                chi_variance_parameter = chi_variance_parameter * pathloss_coefficient;
                std::complex<double> lt_chi = ChiSquare::eval_non_central( -s, non_center_parameter, chi_variance_parameter);
                return lt_chi;
            };

        auto lt_reflected_signal_derivative_s_r = [&](std::complex<double> s, const double y, const double theta){
            double chi_variance_parameter_normalize = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance * direct_link_penalty;
            double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
            double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(test_coverage_hole, y, theta);
            
            double chi_variance_parameter = chi_variance_parameter_normalize * pathloss_coefficient;
            std::complex<double> denominator =  1.0 - 2.0 * chi_variance_parameter * s;
            std::complex<double> lt_chi_D1 = std::exp( s * non_center_parameter * chi_variance_parameter / denominator )/( std::sqrt(denominator) * denominator * denominator);
            lt_chi_D1 = lt_chi_D1 * (1.0 + non_center_parameter - 2.0 * s * chi_variance_parameter)*\
                (chi_variance_parameter * Derivative::partial_argument_BS(test_coverage_hole, BS_density, normalization) + \
                s * chi_variance_parameter_normalize * PathLoss::direct_los(y) * PathLoss::reflected_los_derivative(test_coverage_hole, y, theta) * Derivative::partial_CH_center_BS(BS_density, normalization));
            
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
                PathLoss::reflected_los_derivative(test_coverage_hole, y, theta)*Derivative::partial_CH_center_BS(BS_density, normalization));
        
            return lt_chi_D1; 
        };

        PoissonPP2D_Reduced laplace_transform_interference(lt_interference_at_r, BS_density, test_coverage_hole);
        PoissonPP2D laplace_transform_reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);
        PoissonPP2D_partial laplace_transform_reflected_signals_partial_s_r(lt_reflected_signal_derivative_s_r, RIS_density, cluster_cube);
        PoissonPP2D_partial laplace_transform_reflected_signals_partial_r(lt_reflected_signal_derivative_r, RIS_density, cluster_cube);

        double ergodic_rate;
        double derivative_BS;

        if (RIS_density==0.0){
            std::cout << "No RIS is there " << std::endl;
            ergodic_rate = ErgodicRate::calculate_rayleigh(test_coverage_hole, laplace_transform_interference, noise_power*direct_link_penalty);
            
            derivative_BS = Derivative::get_BS_derivative_rate_no_RIS(test_coverage_hole, BS_density, normalization, \
                noise_power * direct_link_penalty, laplace_transform_interference);
            // std::cout << "ergodic rate: " << ergodic_rate << " derivative_BS " << derivative_BS <<std::endl; 
        } else {
            std::cout << "With RIS " << std::endl;
            ergodic_rate = ErgodicRate::calculate_rayleigh_positive(test_coverage_hole, laplace_transform_interference, noise_power*direct_link_penalty, laplace_transform_reflected_signals);
            derivative_BS = Derivative::get_BS_derivative_rate(test_coverage_hole, BS_density, normalization, \
                noise_power * direct_link_penalty, laplace_transform_interference, laplace_transform_reflected_signals, \
                laplace_transform_reflected_signals_partial_s_r, laplace_transform_reflected_signals_partial_r);
        }


        // double derivative_RIS; 
        double derivative_RIS = Derivative::get_RIS_derivative_rate(test_coverage_hole, noise_power * direct_link_penalty, \
                laplace_transform_interference, laplace_transform_reflected_signals); 


        double cost_of_current_BS = (1.0+RIS_density*(M_PI*(radius_max*radius_max-radius_min*radius_min))/cost_ratio);

        double normalized_gain_BS = derivative_BS / 1e6 /cost_of_current_BS;
        double normailzed_gain_RIS = derivative_RIS / (M_PI*(radius_max*radius_max-radius_min*radius_min)) /(BS_density*1e6) * cost_ratio;

        double RIS_increment = BS_increment * 1e6 * cost_ratio / ( M_PI*( radius_max * radius_max - radius_min * radius_min ) )/(BS_density*1e6);

        BS_gain_vec.push_back(normalized_gain_BS);
        RIS_gain_vec.push_back(normailzed_gain_RIS);
        ergo_rate_vec.push_back(ergodic_rate);



        if (normalized_gain_BS >= normailzed_gain_RIS){
            std::cout << "In this round, we invest in BS!" << std::endl;
            decision_vec.push_back(true);
            BS_density = BS_density + BS_increment/cost_of_current_BS; // Each round, invest two BSs or the same amount
        } else {
            std::cout << "In this round, we invest in RIS!" << std::endl;
            decision_vec.push_back(false);
            RIS_density = RIS_density + RIS_increment;
        }

    }



    json ergo_rate_output(ergo_rate_vec);
    json gain_BS_output(BS_gain_vec);
    json gain_RIS_output(RIS_gain_vec);
    json decision_output(decision_vec);
    json BS_density_output(BS_density_vec);
    json RIS_density_output(RIS_density_vec);

    std::ofstream out_file( "output.json" );
    out_file << "output_rate = " << ergo_rate_output << ";" <<std::endl;
    out_file << "BS_improvement = " << gain_BS_output << ";" <<std::endl;
    out_file << "RIS_improvement = " << gain_RIS_output << ";" <<std::endl;
    out_file << "investment_decision = " << decision_output << ";" <<std::endl;
    out_file << "BS_density = " << BS_density_output << ";" <<std::endl;
    out_file << "RIS_dnesity = " << RIS_density_output << ";" <<std::endl;

    return 0;

}