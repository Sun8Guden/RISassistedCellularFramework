#pragma once
#include <complex>
#include <cmath>
#include "../StochasticGeometry/Parameter.hpp"
#include "../Utils/pathloss.hpp"
#include "../LaplaceTransform/LaplaceTransform.hpp"
#include "../LaplaceTransform/GaussianNoise.hpp"
#include "../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#include "../Integration/boost/math/quadrature/exp_sinh.hpp"
#include "../Integration/GenzMalik/Cube.hpp"
#include "../Integration/GenzMalik/GM2D.hpp"
#include "../../Integration/Cuhre/Integration.hpp"
#include "../../Integration/Cuhre/Region.hpp"

#include <functional>

using namespace boost::math::quadrature;


// Obtaining the derivative w.r.t. BS density for the coverage hole scenrio deserves special attention
// since both the argument and integration region are influenced by the BS density. 

class PoissonPP2D_partial : public LaplaceTransform {
    std::function<std::complex<double>(std::complex<double>, double, double)> mf_laplace_transform_2D;
    CUBE::Cube<double, 2> m_cluster_cube; 
    double m_density_pp; 

    public:
    PoissonPP2D_partial(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D, \
            double density_pp, const CUBE::Cube<double, 2>& cluster_cube): mf_laplace_transform_2D(laplace_function_2D),\
            m_density_pp(density_pp), m_cluster_cube(cluster_cube){}

        void set_functor(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D){
            mf_laplace_transform_2D = laplace_function_2D;
        }

        void set_region(const CUBE::Cube<double, 2>& cluster_cube){
            m_cluster_cube = cluster_cube;
        }
    
    std::complex<double> eval(const std::complex<double>& s) const override {

        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        auto f_integrand = [&](const double y, const double theta){
            return y * ( - mf_laplace_transform_2D(s, y, theta));
        };

        std::complex<double> integral_value = GM::GM2D<double>::integrate(f_integrand, \
            m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return  - integral_value * m_density_pp;
    }

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    std::complex<double> eval_exponent(const std::complex<double>& s) const override {

        std::cout<< "This function is designed to align with the API of the Laplace transform for consistency. Not for use." << std::endl;
        std::complex<double> integral_value {0.0, 0.0};
        return  - integral_value;
    }

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};

namespace Derivative
{

    double partial_CH_center_BS(double BS_density, double normalization_factor){
        return  - normalization_factor / (2.0 * BS_density*std::sqrt(BS_density));
    }

    double partial_argument_BS(double distance, double BS_density, double normalization_factor){
        double return_value = - partial_CH_center_BS(BS_density, normalization_factor) / (Parameter::transmission_power *\
            PathLoss::direct_nlos(distance) * PathLoss::direct_nlos(distance)) * PathLoss::partial_direct_nlos(distance);
        return return_value;
    }

    double partial_noise_CH(const double noise, const double threshold){
        double lt_noise = - noise * threshold / Parameter::antenna_gain; 
        return lt_noise;
    }


    std::complex<double> partial_RIS_CH(std::complex<double> s, \
        const LaplaceTransform& reflected_signals_partial, const double distance){
        std::complex<double> lt_reflected_signal_derivative = reflected_signals_partial.eval( s / PathLoss::direct_nlos(distance) );
        return lt_reflected_signal_derivative;
    }


    std::complex<double> partial_BS_exponent_CH(std::complex<double> s, const double distance, const double threshold, \
        const double BS_density, const double normalization_factor){
        std::complex<double> lower_boundary =  2.0 * M_PI * distance *  threshold * s \
        /(1.0 + s * threshold ) * partial_CH_center_BS(BS_density, normalization_factor);
        double estimated_error{0.0};
        auto f_partial_BS_exponent = [&](const double x){
            double attenuation = threshold * PathLoss::direct_nlos(x);
            std::complex return_value = x * attenuation * Parameter::transmission_power / (1.0 + s * attenuation/ PathLoss::direct_nlos(distance)) / (1.0 + s * attenuation/ PathLoss::direct_nlos(distance));
            return return_value;
        };
        std::complex<double> integral_part = gauss_kronrod<double, 15>::integrate(f_partial_BS_exponent, \
        distance,  Parameter::infinity_large_distance, 10, 1e-6, &estimated_error);
        
        // std::cout<< "lower_boundary" << lower_boundary << " other " << 2.0 * M_PI * integral_part * partial_argument_BS(distance, BS_density, normalization_factor) << "\n";

        return lower_boundary - 2.0 * M_PI * integral_part * partial_argument_BS(distance, BS_density, normalization_factor);
    }

    std::complex<double> partial_BS_exponent_CH_second(std::complex<double> s, const double distance, const double threshold, \
        const double BS_density, const double normalization_factor){
        std::complex<double> lower_boundary =  2.0 * M_PI * distance *  threshold * s \
        /(1.0 + s * threshold ) * partial_CH_center_BS(BS_density, normalization_factor);
        return lower_boundary;
    }

    std::complex<double> partial_BS_CH(std::complex<double> s, const double distance, const double threshold, \
        const double BS_density, const double normalization_factor, \
        const LaplaceTransform& interference){
        
        std::complex<double> lt_interference_exponent = interference.eval_exponent(s * threshold / PathLoss::direct_nlos(distance));
        std::complex<double> lt_interference_derivative = partial_BS_exponent_CH(s, distance, threshold, BS_density, normalization_factor);
        return lt_interference_exponent + BS_density * lt_interference_derivative;
    }

    std::complex<double> partial_BS_CH_second(std::complex<double> s, const double distance, const double threshold, \
        const double BS_density, const double normalization_factor, \
        const LaplaceTransform& interference){
        
        std::complex<double> lt_interference_exponent = interference.eval_exponent(s * threshold/ PathLoss::direct_nlos(distance));
        std::complex<double> lt_interference_derivative = partial_BS_exponent_CH_second(s, distance, threshold, BS_density, normalization_factor);
        return lt_interference_exponent + BS_density * lt_interference_derivative;
    }

    std::complex<double> partial_laplace_CH(std::complex<double> s, \
            const double BS_density, const double normalization_factor, \
            const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals, const LaplaceTransform& reflected_signals_partial_s_r ){
        std::complex<double> lt_value_interference = interference.eval( s * threshold / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_reflected_signal = reflected_signals.eval( s / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);

        std::complex<double> return_value_interference =  lt_value_interference * lt_value_noise * lt_value_reflected_signal *\
            partial_BS_CH(s, distance, threshold, BS_density, normalization_factor, interference);
        std::complex<double> return_value_reflected_signal = lt_value_reflected_signal *  lt_value_interference * lt_value_noise *\
            partial_RIS_CH(s, reflected_signals_partial_s_r, distance);
        std::complex<double> return_value_noise =  lt_value_noise * lt_value_interference * lt_value_reflected_signal * \
            partial_noise_CH(noise, threshold) * partial_argument_BS(distance, BS_density, normalization_factor);

        return  return_value_noise + return_value_interference + return_value_reflected_signal;

    }

    std::complex<double> partial_laplace_CH_no_RIS(std::complex<double> s, \
            const double BS_density, const double normalization_factor, \
            const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference){
        std::complex<double> lt_value_interference = interference.eval( s * threshold / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);

        std::complex<double> return_value_interference =  lt_value_interference * lt_value_noise  * \
            partial_BS_CH(s, distance, threshold, BS_density, normalization_factor, interference);
        std::complex<double> return_value_noise =  lt_value_noise * lt_value_interference * \
            partial_noise_CH(noise, threshold) * partial_argument_BS(distance, BS_density, normalization_factor);

        return  return_value_noise + return_value_interference;

    }

    std::complex<double> partial_laplace_CH_second(std::complex<double> s, \
            const double BS_density, const double normalization_factor, \
            const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals, const LaplaceTransform& reflected_signals_partial_r ){
        std::complex<double> lt_value_interference = interference.eval( s * threshold / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_reflected_signal = reflected_signals.eval( s / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_noise = GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);

        std::complex<double> return_value_interference =  lt_value_interference * lt_value_noise  * lt_value_reflected_signal * \
            partial_BS_CH_second(s, distance, threshold, BS_density, normalization_factor, interference);
        std::complex<double> return_value_reflected_signal = lt_value_reflected_signal *  lt_value_interference * lt_value_noise *\
            partial_RIS_CH(s, reflected_signals_partial_r, distance);
    
        return  return_value_interference + return_value_reflected_signal;

    }




    double get_BS_coverage_no_RIS(const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const double BS_density, const double normalization_factor){
        std::complex<double> lt_value_interference = interference.eval(threshold / PathLoss::direct_nlos(distance));
        std::complex<double> lt_value_noise = GaussianNoise::eval(threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
        std::complex<double> return_value_interference = lt_value_interference * lt_value_noise *\
            partial_BS_CH(std::complex<double>{1.0, 0.0}, distance, threshold, BS_density, normalization_factor, interference);
        std::complex<double> return_value_noise =  lt_value_noise * lt_value_interference * \
            partial_noise_CH(noise, threshold) * partial_argument_BS(distance, BS_density, normalization_factor);
        

        double improvement = return_value_interference.real() + return_value_noise.real();
        return improvement;
    }

    double get_BS_coverage(const double distance, const double threshold, \
            const double BS_density, const double normalization_factor,\
            const double noise, const LaplaceTransform& interference, \
            const LaplaceTransform& reflected_signals, \
            const LaplaceTransform& reflected_signals_partial_s_r, \
            const LaplaceTransform& reflected_signals_partial_r){


        boost::math::quadrature::exp_sinh<double> singular;
        auto f_laplace_eval = [&](const std::complex<double>& s){
            std::complex<double> return_value = partial_laplace_CH(s, BS_density, normalization_factor, distance, threshold, noise, \
                interference, reflected_signals, reflected_signals_partial_s_r);
            return return_value;
        };

        auto f_laplace_second = [&](const std::complex<double>& s){
            std::complex<double> return_value = partial_laplace_CH_second(s, BS_density, normalization_factor, distance, threshold, noise, \
                interference, reflected_signals, reflected_signals_partial_r);
            return return_value;
        };

        auto f_cauchy_principal_value = [&](const double& dummy_z){

            std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, - dummy_z});
            std::complex<double> cpv_part_two = f_laplace_second(std::complex<double>{0.0, - dummy_z});
            std::complex<double> cpv = cpv_part_one - cpv_part_two;
            double coverage_probability_cpv = cpv.imag() / (dummy_z * M_PI);

            return coverage_probability_cpv;
        };
        double termination {Parameter::expected_error};
        double error;
        double L1;
        double cpv_part_value = singular.integrate(f_cauchy_principal_value, termination, &error, &L1);

        std::complex<double> real_part_compute = f_laplace_eval( std::complex<double> {1.0, 0.0});
        double real_part_value = 0.5 * real_part_compute.real();


        return  real_part_value + cpv_part_value;
    }

    double get_BS_derivative_rate_no_RIS(const double distance, \
        const double BS_density, const double normalization_factor,\
        const double noise, const LaplaceTransform& interference){
        
        auto f_laplace_eval = [&](const double threshold){
                std::complex<double> return_value = partial_laplace_CH_no_RIS(std::complex<double>{1.0, 0.0}, BS_density, normalization_factor, distance, threshold, noise, \
                    interference);
                return return_value.real()/(1.0+threshold);
            };

        double estimated_error;
        double derivative_value = gauss_kronrod<double, 15>::integrate(f_laplace_eval, 0.0,  Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error);
        return derivative_value;
    }

    double get_BS_derivative_rate(const double distance, \
        const double BS_density, const double normalization_factor,\
        const double noise, const LaplaceTransform& interference, \
        const LaplaceTransform& reflected_signals, \
        const LaplaceTransform& reflected_signals_partial_s_r, \
        const LaplaceTransform& reflected_signals_partial_r){
            
            std::array<double, 2> center2d{0.0, Parameter::infinity_large_threshold/2.0};
            std::array<double, 2> width2d{Parameter::infinity_cauchy_principal_value,  Parameter::infinity_large_threshold};
            Region<2> region2d = Region<2>(center2d, width2d);
            IntegrationStrategy parallel = IntegrationStrategy::parallel;

            auto f_laplace_eval = [&](const std::complex<double>& s, const double threshold){
                std::complex<double> return_value = partial_laplace_CH(s, BS_density, normalization_factor, distance, threshold, noise, \
                    interference, reflected_signals, reflected_signals_partial_s_r)/(1.0+threshold);
                return return_value;
            };

            auto f_laplace_second = [&](const std::complex<double>& s, const double threshold){
                std::complex<double> return_value = partial_laplace_CH_second(s, BS_density, normalization_factor, distance, threshold, noise, \
                    interference, reflected_signals, reflected_signals_partial_r)/(1.0+threshold);
                return return_value;
            };

            auto f_cpv = [&](const std::array<double, 2>& input2d){
                double z = - std::exp(input2d.at(0));
                std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, z}, input2d.at(1) );
                std::complex<double> cpv_part_two = f_laplace_second(std::complex<double>{0.0, z}, input2d.at(1) );
                std::complex<double> cpv = cpv_part_one - cpv_part_two;
                double coverage_probability_cpv = cpv.imag() / M_PI;
                return coverage_probability_cpv;
            };

            double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
            double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
            int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

            cpv_part_value = Integration<2>::integrate(f_cpv, region2d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 5000, 50000);

            auto f_real = [&](const double threshold){
                std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, threshold);
                double real_part_value = 0.5 * real_part_compute.real();
                return real_part_value;
            };
            real_part_value = gauss_kronrod<double, 15>::integrate(f_real, 0.0,  Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error_part_2);
            

            return cpv_part_value + real_part_value;
        }

} 