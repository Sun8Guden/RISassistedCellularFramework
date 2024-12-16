#pragma once

#include <complex>
#include <iostream>
#include <cmath>
#include <functional>

#include "../../StochasticGeometry/Parameter.hpp"
#include "../../LaplaceTransform/LaplaceTransform.hpp"
#include "../../Integration/Cuhre/Integration.hpp"
#include "../../Integration/Cuhre/Region.hpp"



#define _USE_MATH_DEFINES // for using PI in C++

// Again, this one is a bit special 



namespace ErgodicRateRegion {

        double calculate_rayleigh(const std::function<std::complex<double> (std::complex<double>, double, double)>& f_laplace_eval, const double BS_density){
                
                std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
                std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
                Region<2> region2d = Region<2>(center2d, width2d);
                IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double  real_part_value{ 0.0 };
        double  estimated_error_part_2{ 0.0 };
        int  number_f_eval_part2{ 0 };

                auto f_real = [&](const std::array<double, 2>& input2d ){
                        std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, input2d.at(0), input2d.at(1));
                        double real_part = real_part_compute.real();
                        real_part = real_part * 2.0 * M_PI * BS_density * input2d.at(0) * \
                                std::exp(- M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
                        return real_part;
                };
                real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 50000);

                return real_part_value;
        }

        double calculate_rayleigh_density(const std::function<std::complex<double> (std::complex<double>, double, double)>& f_laplace_eval, const double BS_density){
                
                std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
                std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
                Region<2> region2d = Region<2>(center2d, width2d);
                IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double  real_part_value{ 0.0 };
        double  estimated_error_part_2{ 0.0 };
        int  number_f_eval_part2{ 0 };

                auto f_real = [&](const std::array<double, 2>& input2d ){
                        std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, input2d.at(0), input2d.at(1));
                        double real_part = real_part_compute.real();
                        real_part = real_part * 2.0 * M_PI * input2d.at(0) * (1.0 - M_PI * BS_density * input2d.at(0) * input2d.at(0)) * \
                                std::exp(- M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
                        return real_part;
                }; 
                real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 50000);

                return real_part_value;
        }


    double calculate_rayleigh_positive(const std::function<std::complex<double> (std::complex<double>, double, double)>& f_laplace_eval, const double BS_density){
        std::array<double, 3> center3d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0,  0.0};
        std::array<double, 3> width3d{Parameter::infinity_large_distance-Parameter::guard_distance, Parameter::infinity_large_threshold,  Parameter::infinity_cauchy_principal_value};
        Region<3> region3d = Region<3>(center3d, width3d);

        std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
        std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
        Region<2> region2d = Region<2>(center2d, width2d);
        IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
        double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
        int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

        auto f_cpv = [&](const std::array<double, 3>& input3d) {
                double z = - std::exp(input3d.at(2));
                std::complex<double> cpv_part_one = f_laplace_eval( std::complex<double>{1.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv_part_two = f_laplace_eval( std::complex<double>{0.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv =  cpv_part_one - cpv_part_two;  
                double cpv_val = cpv.imag();

                double ret_val =  cpv_val * 2.0 * BS_density * input3d.at(0) * \
                        std::exp(- M_PI * BS_density * input3d.at(0) * input3d.at(0))/(1.0+input3d.at(1));

                return ret_val;
        };

        cpv_part_value = Integration<3>::integrate(f_cpv, region3d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 5000, 1000000);

        std::cout << " evluation " << number_f_eval_part1  << " error " << estimated_error_part_1 << std::endl;

        auto f_real = [&](const std::array<double, 2>& input2d ){
            std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, input2d.at(0), input2d.at(1));
            double real_part =  0.5 * ( 1.0 + real_part_compute.real());
            real_part = real_part * 2.0 * M_PI * BS_density * input2d.at(0) * \
                        std::exp(- M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
            return real_part;
        };

        real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 50000);
        std::cout << " evluation 2 " << number_f_eval_part2  << " error " << estimated_error_part_2 << std::endl;

        // std::cout << "CPV " << cpv_part_value << "real " << real_part_value <<"\n";
        double estimated_ergodic_rate = real_part_value + cpv_part_value;

        

        return estimated_ergodic_rate;
    }

    double calculate_from_probability(const std::function<double (double, double)>& f_laplace_eval, const double BS_density){
                std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
                std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
                Region<2> region2d = Region<2>(center2d, width2d);
                IntegrationStrategy parallel = IntegrationStrategy::parallel;
                double  real_part_value{ 0.0 };
                double  estimated_error_part_2{ 0.0 };
                int  number_f_eval_part2{ 0 };

                auto f_real = [&](const std::array<double, 2>& input2d ){
                        std::complex<double> real_part_compute = f_laplace_eval(input2d.at(0), input2d.at(1));
                        double real_part = real_part_compute.real();
                        real_part = real_part * 2.0 * M_PI * BS_density * input2d.at(0) * \
                                std::exp(- M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
                        return real_part;
                };
                real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 50000);

                return real_part_value;
    }

    double calculate_rayleigh_positive_derivative(const std::function<std::complex<double> (std::complex<double>, double, double)>& f_laplace_eval, const double BS_density){
        std::array<double, 3> center3d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0,  0.0};
        std::array<double, 3> width3d{Parameter::infinity_large_distance-Parameter::guard_distance, Parameter::infinity_large_threshold,  Parameter::infinity_cauchy_principal_value};
        Region<3> region3d = Region<3>(center3d, width3d);

        std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
        std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
        Region<2> region2d = Region<2>(center2d, width2d);
        IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
        double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
        int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

        auto f_cpv = [&](const std::array<double, 3>& input3d) {
                double z = - std::exp(input3d.at(2));
                std::complex<double> cpv_part_one = f_laplace_eval( std::complex<double>{1.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv_part_two = f_laplace_eval( std::complex<double>{0.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv =  cpv_part_one - cpv_part_two;  
                double cpv_val = cpv.imag();

                double ret_val =  cpv_val * 2.0 * BS_density * input3d.at(0) * \
                        std::exp(- M_PI * BS_density * input3d.at(0) * input3d.at(0))/(1.0+input3d.at(1));

                return ret_val;
        };

        cpv_part_value = Integration<3>::integrate(f_cpv, region3d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 5000, 1000000);

        std::cout << " evluation " << number_f_eval_part1  << " error " << estimated_error_part_1 << std::endl;

        auto f_real = [&](const std::array<double, 2>& input2d ){
            std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, input2d.at(0), input2d.at(1));
            double real_part =  0.5 * ( real_part_compute.real() );
            real_part = real_part * 2.0 * M_PI * BS_density * input2d.at(0) * \
                        std::exp(- M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
            return real_part;
        };

        real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 40000);
        std::cout << " evluation 2 " << number_f_eval_part2  << " error " << estimated_error_part_2 << std::endl;

        // std::cout << "CPV " << cpv_part_value << "real " << real_part_value <<"\n";
        double estimated_ergodic_rate = real_part_value + cpv_part_value;

        

        return estimated_ergodic_rate;

    }


    double calculate_rayleigh_positive_density_derivative(const std::function<std::complex<double> (std::complex<double>, double, double)>& f_laplace_eval, const double BS_density){
        std::array<double, 3> center3d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0,  0.0};
        std::array<double, 3> width3d{Parameter::infinity_large_distance-Parameter::guard_distance, Parameter::infinity_large_threshold,  Parameter::infinity_cauchy_principal_value};
        Region<3> region3d = Region<3>(center3d, width3d);

        std::array<double, 2> center2d{(Parameter::infinity_large_distance+Parameter::guard_distance)/2.0, Parameter::infinity_large_threshold/2.0};
        std::array<double, 2> width2d{Parameter::infinity_large_distance-Parameter::guard_distance,  Parameter::infinity_large_threshold};
        Region<2> region2d = Region<2>(center2d, width2d);
        IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
        double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
        int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

        auto f_cpv = [&](const std::array<double, 3>& input3d) {
                double z = - std::exp(input3d.at(2));
                std::complex<double> cpv_part_one = f_laplace_eval( std::complex<double>{1.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv_part_two = f_laplace_eval( std::complex<double>{0.0, z}, input3d.at(0), input3d.at(1) );
                std::complex<double> cpv =  cpv_part_one - cpv_part_two;  
                double cpv_val = cpv.imag();

                double ret_val =  cpv_val * ( 2.0 * input3d.at(0) - 2.0 * M_PI * input3d.at(0) * input3d.at(0) * input3d.at(0) * BS_density) * \
                        std::exp(- M_PI * BS_density * input3d.at(0) * input3d.at(0)) / (1.0+input3d.at(1));

                return ret_val;
        };

        cpv_part_value = Integration<3>::integrate(f_cpv, region3d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 5000, 1000000);

        std::cout << " evluation " << number_f_eval_part1  << " error " << estimated_error_part_1 << std::endl;

        auto f_real = [&](const std::array<double, 2>& input2d ){
            std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, input2d.at(0), input2d.at(1));
            double real_part =  0.5 * (1.0 + real_part_compute.real() );
            real_part = real_part * M_PI * ( 2.0 * input2d.at(0) - 2.0 * M_PI * input2d.at(0) * input2d.at(0) * input2d.at(0) * BS_density) * \
                        std::exp( - M_PI * BS_density * input2d.at(0) * input2d.at(0))/(1.0+input2d.at(1));
            return real_part;
        };

        real_part_value = Integration<2>::integrate(f_real, region2d, estimated_error_part_2, \
                number_f_eval_part2, parallel, 1e-6, 500, 40000);
        std::cout << " evluation 2 " << number_f_eval_part2  << " error " << estimated_error_part_2 << std::endl;

        // std::cout << "CPV " << cpv_part_value << "real " << real_part_value <<"\n";
        double estimated_ergodic_rate = real_part_value + cpv_part_value;

        return estimated_ergodic_rate;

    }


}