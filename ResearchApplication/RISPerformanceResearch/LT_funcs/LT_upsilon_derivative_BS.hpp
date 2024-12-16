#pragma once
#include "LaplaceTransform.hpp"
#include "LT_interference.hpp"
#include "LT_signal_chi.hpp"
#include "LT_signal_gamma.hpp"
#include "library/GenzMalik/Cube.hpp"
#include "StochasticGeometryArchiv/Param.hpp"

class LT_upsilon_derivative_BS {

    LT_interference interference;
    LT_signal_chi signal_chi;
    Param param;
    // LT_signal_gamma signal_gamma;

    public:
    LT_upsilon_derivative_BS(const Param& param_, CUBE::Cube<double, 2>& cluster_cube_): interference(param_), signal_chi(param_, cluster_cube_), param(param_){}
    ~LT_upsilon_derivative_BS(){}

    std::array<std::complex<double>, 2> lt_eval_array(const std::complex<double>& s, const double& threshold, const double& distance) {
        std::array<std::complex<double>, 2> eval_result;

        std::complex<double> interference_eval =  interference.lt_eval(s, threshold, distance);
        std::complex<double> signal_eval {1.0, 0.0}; 
        if (param.RISden != 0.0){
            signal_eval = signal_chi.lt_eval(s, threshold, distance);
        }
        std::complex<double> noise_term = s / PL_d(distance) / param.Antenna_gain * param.noise_level * threshold;
        std::complex<double> exp_term = std::exp(- param.BSden * interference_eval - param.RISden * signal_eval  - noise_term);
        eval_result.at(0) = exp_term * ( - interference_eval );
        eval_result.at(1) = exp_term;
        return eval_result;
 
    }

};