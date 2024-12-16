#pragma once
#include <complex>
#include <iostream>
#include <cmath>
#include <boost/math/quadrature/exp_sinh.hpp>
using namespace boost::math::quadrature;
exp_sinh<double> singular;

template <typename LT>
double separate_positive(LT lt_interference_signal, std::complex<double> s_parameter, double distance_arg, double threshold_arg)
{

    auto f_cpv = [&](const double& dummy_z){
        std::complex<double> lt_vec_part_1_one = lt_interference_signal(s_parameter - std::complex<double>(0.0, 1.0) * dummy_z, distance_arg, threshold_arg);
        std::complex<double> lt_vec_part_1_zero = lt_interference_signal( - std::complex<double>(0.0, 1.0) * dummy_z, distance_arg, threshold_arg);

        std::complex<double> cpv = (lt_vec_part_1_one - lt_vec_part_1_zero);
        double cpv_val = cpv.imag() / (dummy_z * M_PI);
        // std::cout << "VALUE Z "<< dummy_z << " VALUE ONE " << lt_vec_part_1_one << std::endl;
        return cpv_val;
    };

    double termination {1e-6};
    double error;   double L1;
    double cpv = singular.integrate(f_cpv, termination, &error, &L1);
    double value_part_1 = cpv;

    std::complex<double> lt_part2 = 1.0 + lt_interference_signal(s_parameter, distance_arg, threshold_arg);
    double value_part_2 = 0.5 * lt_part2.real();
    // std::cout << "value_part_1 " << value_part_1 << "value_part_2 " << value_part_2<< std::endl;

    return (value_part_1 + value_part_2);
}
