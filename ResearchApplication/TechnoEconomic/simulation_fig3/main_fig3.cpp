#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>

#include "../RateRegion.hpp"

#include "../../../LaplaceTransform/LaplaceTransform.hpp"
#include "../../../LaplaceTransform/PoissonPP.hpp"
#include "../../../LaplaceTransform/Exponential.hpp"
#include "../../../LaplaceTransform/GaussianNoise.hpp"
#include "../../../LaplaceTransform/ChiSquare.hpp"
#include "../../../Utils/pathloss.hpp"


#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"
#include "../../../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#define _USE_MATH_DEFINES // for using PI in C++

#include "../../../Utils/json.hpp"
using json = nlohmann::json;

int main(){

    // Each round, one has the resource to deploy one additional BS per km^2, we show the route of deployment strategy by varying the cost ratio!
    double cost_ratio = 5.0;
    double BS_increment = 2e-6;

    double noise_power = 1e-13; // Watt 
    double RIS_elements = 600.0;
    double radius_min=20.0;
    double radius_max=30.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    double BS_density = 1e-6; // Initial BS density.
    double RIS_density = {0.0/ (M_PI*(radius_max*radius_max-radius_min*radius_min))};



    auto lt_interference_at_r = [&](std::complex<double> s, double r){
        double interference_power = PathLoss::direct_nlos(r); 
        std::complex<double> return_value = Exponential::eval(s, interference_power);
        return return_value;
    };

    auto lt_reflected_signal_at_x_y_theta = [&](std::complex<double> s, const double x, const double y, const double theta){
        double chi_variance_parameter = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance;
        double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean /Parameter::zeta_variance;
        double pathloss_coefficient = PathLoss::direct_los(y) * PathLoss::reflected_los(x, y, theta);
        chi_variance_parameter = chi_variance_parameter * pathloss_coefficient;
        std::complex<double> lt_chi = ChiSquare::eval_non_central( - s, non_center_parameter, chi_variance_parameter);
        return lt_chi;
    };

    std::vector<double> ergo_rate_vec;
    std::vector<double> BS_gain_vec;
    std::vector<double> RIS_gain_vec;
    std::vector<bool>   decision_vec;
    std::vector<double> BS_density_vec;
    std::vector<double> RIS_density_vec;


    for (int investment_round = 1; investment_round < 20; investment_round++){


        std::cout << "Investment Round: " << investment_round << " Round!\n";
        std::cout << "BS density is " << BS_density << ", RIS density is " <<  RIS_density <<  std::endl;
        BS_density_vec.push_back(BS_density);
        RIS_density_vec.push_back(RIS_density);

        auto lt_ergodic_rate = [&](std::complex<double> s, const double distance, const double threshold){
            double estimated_error_BS {0.0};
            auto f_integrand_interference = [&](const double x){
                return x * (1.0 - lt_interference_at_r( s * threshold  / PathLoss::direct_nlos(distance), x));
            };
            std::complex<double> lt_interference_exponent = gauss_kronrod<double, 15>::integrate(f_integrand_interference, \
            distance,  Parameter::infinity_large_distance, 10, 1e-6, &estimated_error_BS);
            std::complex<double> lt_interference = std::exp( - 2.0 * M_PI * lt_interference_exponent * BS_density); 


            double estimated_error_RIS{ 0.0 };
            unsigned int number_evaluation{0};
            auto f_integrand_RIS = [&](const double y, const double theta){
                return y * (1.0 - lt_reflected_signal_at_x_y_theta(s / PathLoss::direct_nlos(distance), distance, y, theta));
            };
            std::complex<double> lt_signal_exponent = GM::GM2D<double>::integrate(f_integrand_RIS, \
                cluster_cube, 1e-6, estimated_error_RIS, 5000, number_evaluation);
            std::complex<double> lt_signal = std::exp( - lt_signal_exponent * RIS_density );

            std::complex<double> lt_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise_power);

            std::complex<double> integrand =  lt_interference * lt_noise * lt_signal; 

            return integrand;

        };

        // Laplace transform of derivatives: 
        auto lt_ergodic_rate_derivative_BS = [&](std::complex<double> s, const double distance, const double threshold){
            double estimated_error_BS {0.0};
            auto f_integrand_interference = [&](const double x){
                return x * (1.0 - lt_interference_at_r( s * threshold  / PathLoss::direct_nlos(distance), x));
            };
            std::complex<double> lt_interference_exponent = gauss_kronrod<double, 15>::integrate(f_integrand_interference, \
            distance,  Parameter::infinity_large_distance, 10, 1e-6, &estimated_error_BS);
            std::complex<double> lt_interference = std::exp( - 2.0 * M_PI * lt_interference_exponent * BS_density); 


            double estimated_error_RIS{ 0.0 };
            unsigned int number_evaluation{0};
            auto f_integrand_RIS = [&](const double y, const double theta){
                return y * (1.0 - lt_reflected_signal_at_x_y_theta(s / PathLoss::direct_nlos(distance), distance, y, theta));
            };
            std::complex<double> lt_signal_exponent = GM::GM2D<double>::integrate(f_integrand_RIS, \
                cluster_cube, 1e-6, estimated_error_RIS, 5000, number_evaluation);
            std::complex<double> lt_signal = std::exp( - lt_signal_exponent * RIS_density );

            std::complex<double> lt_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise_power);

            std::complex<double> integrand =  lt_interference * lt_noise * lt_signal * (- 2.0 * M_PI * lt_interference_exponent); 

            return integrand;

        };

        auto lt_ergodic_rate_derivative_RIS = [&](std::complex<double> s, const double distance, const double threshold){
            double estimated_error_BS {0.0};
            auto f_integrand_interference = [&](const double x){
                return x * (1.0 - lt_interference_at_r( s * threshold  / PathLoss::direct_nlos(distance), x));
            };
            std::complex<double> lt_interference_exponent = gauss_kronrod<double, 15>::integrate(f_integrand_interference, \
            distance,  Parameter::infinity_large_distance, 10, 1e-6, &estimated_error_BS);
            std::complex<double> lt_interference = std::exp( - 2.0 * M_PI * lt_interference_exponent * BS_density); 


            double estimated_error_RIS{ 0.0 };
            unsigned int number_evaluation{0};
            auto f_integrand_RIS = [&](const double y, const double theta){
                return y * (1.0 - lt_reflected_signal_at_x_y_theta(s / PathLoss::direct_nlos(distance), distance, y, theta));
            };
            std::complex<double> lt_signal_exponent = GM::GM2D<double>::integrate(f_integrand_RIS, \
                cluster_cube, 1e-6, estimated_error_RIS, 5000, number_evaluation);
            std::complex<double> lt_signal = std::exp( - lt_signal_exponent * RIS_density );
            std::complex<double> lt_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise_power);
            std::complex<double> integrand =  lt_interference * lt_noise * lt_signal * ( - lt_signal_exponent ); 
            return integrand;
        };
                    
        
        double derivative_BS;
        double ergodic_rate;
        if (RIS_density == 0.0){
            derivative_BS = ErgodicRateRegion::calculate_rayleigh(lt_ergodic_rate_derivative_BS, BS_density) + \
                ErgodicRateRegion::calculate_rayleigh_density(lt_ergodic_rate, BS_density);

            ergodic_rate = ErgodicRateRegion::calculate_rayleigh(lt_ergodic_rate, BS_density);
        } else {
            derivative_BS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_BS, BS_density) +\
                ErgodicRateRegion::calculate_rayleigh_positive_density_derivative(lt_ergodic_rate, BS_density);
            ergodic_rate = ErgodicRateRegion::calculate_rayleigh_positive(lt_ergodic_rate, BS_density);
        }
        
        // derivative_BS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_BS, BS_density) +\
        //             ErgodicRateRegion::calculate_rayleigh_positive_density_derivative(lt_ergodic_rate, BS_density);

        double derivative_RIS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_RIS, BS_density);



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
    out_file << "investment decision = " << decision_output << ";" <<std::endl;
    out_file << "BS_density = " << BS_density_output << ";" <<std::endl;
    out_file << "RIS_dnesity = " << RIS_density_output << ";" <<std::endl;


return 0;


}