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

    auto start{std::chrono::steady_clock::now()};


    double noise_power = 1e-13; // Watt 
    double RIS_elements = 200.0;
    double radius_min=20.0;
    double radius_max=30.0;

    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);
    double RIS_density = {0.0/ (M_PI*(radius_max*radius_max-radius_min*radius_min))};


    // Laplace transform preparation

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


    std::vector<double> rate_gain_BS;
    std::vector<double> rate_gain_RIS;



    // Begin the computation!!

    for (int i=1; i<20; i++){


        std::cout << "Running " << i << "Round!\n";
        double BS_density = (double(i) * 2.0 + 1.0) * 1e-6;

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
        if (RIS_density == 0.0){
            derivative_BS = ErgodicRateRegion::calculate_rayleigh(lt_ergodic_rate_derivative_BS, BS_density) + \
                ErgodicRateRegion::calculate_rayleigh_density(lt_ergodic_rate, BS_density);
        }else{
            derivative_BS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_BS, BS_density) +\
                ErgodicRateRegion::calculate_rayleigh_positive_density_derivative(lt_ergodic_rate, BS_density);
        }
        
        // derivative_BS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_BS, BS_density) +\
        //             ErgodicRateRegion::calculate_rayleigh_positive_density_derivative(lt_ergodic_rate, BS_density);



        double derivative_RIS = ErgodicRateRegion::calculate_rayleigh_positive_derivative(lt_ergodic_rate_derivative_RIS, BS_density);

        double normalized_gain_BS = derivative_BS * 1e-6;
        double normailzed_gain_RIS = derivative_RIS / (M_PI*(radius_max*radius_max-radius_min*radius_min)) / (BS_density * 1e6);

        rate_gain_BS.push_back(normalized_gain_BS);
        rate_gain_RIS.push_back(normailzed_gain_RIS);

    }


    json gain_BS_output(rate_gain_BS);
    json gain_RIS_output(rate_gain_RIS);

    std::ofstream out_file( "output.json" );
    out_file << "BS improvement is = " << gain_BS_output << ";" <<std::endl;
    out_file << "RIS improvement is  = " << gain_RIS_output << ";" <<std::endl;

    return 0;

}