#pragma once
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <complex>
#include "Param.hpp"
#include "library/boost_1_83_0/boost/math/quadrature/exp_sinh.hpp"
#include "LaplacetransformArchiv/LaplaceTransform.hpp"
using namespace boost::math::quadrature;

class Cov_Prob {
    exp_sinh<double> singular;
    std::complex<double> imag_unit;
    Param param;
    // The param is here for switching to the corresponding parameters sets. 

    double cov_prob_BF_1(LaplaceTransform& lt_trans, const double& s_BF, const double& threshold, const double& distance_BS){
        std::complex<double> s(s_BF, 0);
        auto f_cpv = [&](const double& dummy_z){
            // Real z = std::exp(z_orig);
            std::complex<double> lt_vec_part_1_one = lt_trans.lt_eval(s - imag_unit * dummy_z, threshold, distance_BS);
            std::complex<double> lt_vec_part_1_zero = lt_trans.lt_eval( - imag_unit * dummy_z, threshold, distance_BS);
            std::complex<double> cpv = (lt_vec_part_1_one - lt_vec_part_1_zero);
            // std::cout << "cpv value is " << cpv << std::endl;
            double cpv_val = cpv.imag() / (dummy_z * M_PI);
            return cpv_val;
        };
        double termination {1e-6};
        double error;
        double L1;
        double cpv = singular.integrate(f_cpv, termination, &error, &L1);
        double part_1_val = cpv;

        std::complex<double> lt_vec_part_2 = lt_trans.lt_eval(s_BF, threshold, distance_BS);
        std::complex<double> lt_part2 = 1.0 + lt_vec_part_2;
        double part_2_val = 0.5 * lt_part2.real();
        return (part_1_val + part_2_val);
    }

public:
        Cov_Prob(const Param& param_): param(param_), imag_unit(0.0, 1.0){}
        ~Cov_Prob(){}

        double calc(LaplaceTransform& lt_trans, const double& threshold, const double& distance_BS){
            double cov_prob {0.0};
            cov_prob = cov_prob_BF_1(lt_trans, 1.0, threshold, distance_BS);
            return cov_prob;
        }

        double calc(LaplaceTransform& lt_trans, const double& s_BF,  const double& threshold, const double& distance_BS){
            double cov_prob {0.0};
            cov_prob = cov_prob_BF_1(lt_trans, s_BF, threshold, distance_BS);
            return cov_prob;
        }

};