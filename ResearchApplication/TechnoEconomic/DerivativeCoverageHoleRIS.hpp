#pragma once
#include <complex>
#include "../../LaplaceTransform/LaplaceTransform.hpp"
#include "../../LaplaceTransform/GaussianNoise.hpp"
#include "../../Utils/pathloss.hpp"
#include "../../StochasticGeometry/Parameter.hpp"
#include "../../Integration/Cuhre/Integration.hpp"
#include "../../Integration/Cuhre/Region.hpp"
#include "../Integration/boost/math/quadrature/gauss_kronrod.hpp"


namespace Derivative {


    double get_RIS_derivative_wo_RIS(const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals){
        double coverage_probabilility = interference.eval( threshold  / PathLoss::direct_nlos(distance) )*\
            GaussianNoise::eval(threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise);
        coverage_probabilility = coverage_probabilility * reflected_signals.eval_exponent(1.0/PathLoss::direct_nlos(distance));
        return coverage_probabilility;
    }

    double get_RIS_derivative(const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals){

        boost::math::quadrature::exp_sinh<double> singular;
        auto f_laplace_eval = [&](const std::complex<double>& s){
            std::complex<double> coverage_probabilility = interference.eval( s * threshold/PathLoss::direct_nlos(distance) )*\
            reflected_signals.eval( s / PathLoss::direct_nlos(distance) )*\
            GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
            coverage_probabilility = coverage_probabilility * reflected_signals.eval_exponent(s /PathLoss::direct_nlos(distance));
            return coverage_probabilility;
        };

        auto f_cauchy_principal_value = [&](const double& dummy_z){
            std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, -dummy_z});
            std::complex<double> cpv_part_two = f_laplace_eval(std::complex<double>{0.0, -dummy_z});
            std::complex<double> cpv = cpv_part_one - cpv_part_two;
            double coverage_probability_cpv = cpv.imag() / (dummy_z * M_PI);
            return coverage_probability_cpv;
        };
        double termination {Parameter::expected_error};
        double error;
        double L1;
        double cpv_part_value = singular.integrate(f_cauchy_principal_value, termination, &error, &L1);

        std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0});
        double real_part_value = 0.5 * real_part_compute.real();

        return cpv_part_value + real_part_value;
    }

    double get_RIS_derivative_rate(const double distance, const double noise, \
        const LaplaceTransform& interference, const LaplaceTransform& reflected_signals){

            std::array<double, 2> center2d{0.0, Parameter::infinity_large_threshold/2.0};
            std::array<double, 2> width2d{Parameter::infinity_cauchy_principal_value,  Parameter::infinity_large_threshold};
            Region<2> region2d = Region<2>(center2d, width2d);
            IntegrationStrategy parallel = IntegrationStrategy::parallel;

            auto f_laplace_eval = [&](const std::complex<double>& s, const double threshold){
                std::complex<double> coverage_probabilility = interference.eval( s * threshold/PathLoss::direct_nlos(distance) )*\
                reflected_signals.eval( s / PathLoss::direct_nlos(distance) )*\
                GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
                coverage_probabilility = coverage_probabilility * reflected_signals.eval_exponent(s/PathLoss::direct_nlos(distance));
                return coverage_probabilility/(1.0+threshold);
            };

            auto f_cpv = [&](const std::array<double, 2>& input2d){
                double z = - std::exp(input2d.at(0));
                std::complex<double> cpv_part_one = f_laplace_eval(std::complex<double>{1.0, z}, input2d.at(1) );
                std::complex<double> cpv_part_two = f_laplace_eval(std::complex<double>{0.0, z}, input2d.at(1) );
                std::complex<double> cpv = cpv_part_one - cpv_part_two;
                double coverage_probability_cpv = cpv.imag() / M_PI;
                return coverage_probability_cpv;
            };
            
            double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
            double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
            int number_f_eval_part1{ 0 }, number_f_eval_part2{ 0 };

            boost::math::quadrature::exp_sinh<double> singular;
            cpv_part_value = Integration<2>::integrate(f_cpv, region2d, estimated_error_part_1, \
                number_f_eval_part1, parallel, 1e-6, 5000, 50000);

            auto f_real = [&](const double threshold){
                std::complex<double> real_part_compute = f_laplace_eval(std::complex<double>{1.0, 0.0}, threshold);
                double real_part_value = 0.5 * real_part_compute.real();
                return real_part_value;
            };

            real_part_value = gauss_kronrod<double, 15>::integrate(f_real, 0.0,  Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error_part_2);

        // std::cout << "CPV part is: " << cpv_part_value << "Real value part: " << real_part_value << std::endl;
        return cpv_part_value + real_part_value; 
    }

}