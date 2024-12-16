#pragma once
#include <complex>
#include <iostream>
#include <cmath>
#include "Parameter.hpp"
#include "CoverageProbability.hpp"
#include "../LaplaceTransform/LaplaceTransform.hpp"
#include "../Integration/Cuhre/Integration.hpp"
#include "../Integration/Cuhre/Region.hpp"
#include "../Integration/GenzMalik/Cube.hpp"
#include "../Integration/GenzMalik/GM2D.hpp"
#include <../Integration/boost/math/quadrature/gauss_kronrod.hpp>
#define _USE_MATH_DEFINES // for using PI in C++
using namespace boost::math::quadrature;

// We will use multiple integral here.

class ErgodicRate {

    public:
    
    static double calculate_rayleigh(const double distance, const LaplaceTransform& interference, const double& noise){
        auto coverage_probability = [&](double threshold){
            double coverage_val =  CoverageProbability::calculate_rayleigh(distance, threshold, interference, noise);
            return coverage_val /(1.0 + threshold);
        };
        double estimated_error;
        double estimated_ergodic_rate = gauss_kronrod<double, 15>::integrate(coverage_probability, 0.0, \
                Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error);
        return estimated_ergodic_rate;
    }

    static double calculate_rayleigh_positive(const double distance, const LaplaceTransform& interference, \
            const double noise, const LaplaceTransform& reflected_signals){

        auto f_laplace_eval = [&](const std::complex<double>& s, const double threshold){
            std::complex<double> coverage_probabilility = interference.eval( s * threshold/PathLoss::direct_nlos(distance) ) *\
            reflected_signals.eval( s / PathLoss::direct_nlos(distance)) * \
            GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
            return coverage_probabilility / (1.0 + threshold);
        };

        std::array<double, 2> center2d{Parameter::infinity_large_threshold/2.0,  0.0};
        std::array<double, 2> width2d{Parameter::infinity_large_threshold,  Parameter::infinity_cauchy_principal_value};
        Region<2> cur_region2d = Region<2>(center2d, width2d);
        IntegrationStrategy parallel = IntegrationStrategy::parallel;

        double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
        double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
        int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

        auto f_cpv = [&](const std::array<double, 2>& input2d) {
            double z = - std::exp(input2d.at(1));
            std::complex<double> cpv_part_one = f_laplace_eval( std::complex<double>{1.0, z}, input2d.at(0) );
            std::complex<double> cpv_part_two = f_laplace_eval( std::complex<double>{0.0, z}, input2d.at(0) );
            std::complex<double> cpv =  cpv_part_one - cpv_part_two;  
            double cpv_val = cpv.imag();
            double ret_val = cpv_val / M_PI;
            return ret_val;
        };
        cpv_part_value = Integration<2>::integrate(f_cpv, cur_region2d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 500, 50000);
        
        auto f_real = [&](const double threshold ){
            std::complex<double> real_part_compute =  f_laplace_eval(std::complex<double>{1.0, 0.0}, threshold);
            double real_part = 0.5 * (1.0/(1.0+threshold) + real_part_compute.real());
            return real_part;

        };
        real_part_value =  gauss_kronrod<double, 15>::integrate(f_real, 0.0, Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error_part_2); 

        // std::cout << "real_part_value " << real_part_value << "cpv_part_value " <<cpv_part_value << std::endl;
        double estimated_ergodic_rate = real_part_value + cpv_part_value;
        return estimated_ergodic_rate;
    }
};
