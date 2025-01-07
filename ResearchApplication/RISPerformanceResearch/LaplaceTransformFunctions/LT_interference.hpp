#pragma once
#include <cmath>
#include <iostream>
#include <complex>
#include "../../../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#include "LaplaceTransform.hpp"
#include "LT_utils.hpp"
#include "../StochasticGeometryHelper/Param.hpp"

class LT_interference : public LaplaceTransform {

    Param param;

    public:
    LT_interference(const Param & param_): param(param_){}
    ~LT_interference(){}

    std::complex<double> lt_eval(const std::complex<double>& s, const double& threshold, const double& distance) override {
        auto f_assist = [&](const double& x){
            std::complex<double> ret_val  = x * (1.0 - 1.0 / (1.0 + s * threshold * PL_d(x) / PL_d(distance)));
            return ret_val;
        };
        double estimated_error_f_assist;
        std::complex<double> out_integral = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f_assist, distance, 5000.0, 10, 1e-6, &estimated_error_f_assist); 
        std::complex<double> lt_interference = 2.0 * M_PI * out_integral;
        return lt_interference;
    }
};