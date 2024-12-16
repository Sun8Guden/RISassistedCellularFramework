#pragma once
#define _USE_MATH_DEFINES // for using the constant PI in C++
#include <cmath>
#include <iostream>
#include <complex>
#include "Parameter.hpp"
#include "../Integration/boost/math/quadrature/exp_sinh.hpp"
#include "../LaplaceTransform/LaplaceTransform.hpp"
#include "../LaplaceTransform/PoissonPP.hpp"
#include "../LaplaceTransform/GaussianNoise.hpp"
#include "../Utils/pathloss.hpp"


using namespace boost::math::quadrature;

// We want to provide a general methodology for computing coverage probability for multiple common direct links,
// with or without the reflected signals 
namespace CoverageProbability {


     
    double calculate_rayleigh(const double distance, const double threshold, \
            const LaplaceTransform& interference, const double noise){
        
        double coverage_probabilility = interference.eval( threshold  / PathLoss::direct_nlos(distance) )*\
            GaussianNoise::eval(threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise);
        //
        double a =  threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain * Parameter::transmission_power);
        // std::cout << "Function called at " << distance <<  " with noise " << a << std::endl;

        return coverage_probabilility;

    }

    double calculate_rayleigh_positive(const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals)
    {
        boost::math::quadrature::exp_sinh<double> singular;

        auto f_laplace_eval = [&](const std::complex<double>& s){
            std::complex<double> coverage_probabilility = interference.eval( s * threshold/PathLoss::direct_nlos(distance) )*\
            reflected_signals.eval( s / PathLoss::direct_nlos(distance) ) *\
            GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
            //
            return coverage_probabilility;
        };

        auto f_cauchy_principal_value = [&](const double& dummy_z){
            std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, - dummy_z});
            std::complex<double> cpv_part_two = f_laplace_eval(std::complex<double>{0.0, - dummy_z});
            std::complex<double> cpv = cpv_part_one - cpv_part_two;
            double coverage_probability_cpv = cpv.imag() / (dummy_z * M_PI);
            return coverage_probability_cpv;
        };
        double termination {Parameter::expected_error};
        double error;
        double L1;
        double cpv_part_value = singular.integrate(f_cauchy_principal_value, termination, &error, &L1);

        std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0});
        double real_part_value = 0.5 * (1.0 + real_part_compute.real() );
        double coverage_probabilility = cpv_part_value + real_part_value;
        return coverage_probabilility;
    }


    double calculate_laplace_positive(std::function<std::complex<double>(std::complex<double>, double, double)> f_laplace_eval, const double distance, const double threshold){
        boost::math::quadrature::exp_sinh<double> singular;

        auto f_cauchy_principal_value = [&](const double& dummy_z){
            std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, - dummy_z}, distance, threshold);
            std::complex<double> cpv_part_two = f_laplace_eval(std::complex<double>{0.0, - dummy_z}, distance, threshold);
            std::complex<double> cpv = cpv_part_one - cpv_part_two;
            double coverage_probability_cpv = cpv.imag() / (dummy_z * M_PI);
            return coverage_probability_cpv;
        };
        double termination {Parameter::expected_error};
        double error;
        double L1;
        double cpv_part_value = singular.integrate(f_cauchy_principal_value, termination, &error, &L1);

        std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, distance, threshold);
        double real_part_value = 0.5 * (1.0 + real_part_compute.real() );
        double coverage_probabilility = cpv_part_value + real_part_value;
        return coverage_probabilility;

    }
}