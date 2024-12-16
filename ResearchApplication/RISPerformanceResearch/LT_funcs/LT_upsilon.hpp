#pragma once
#include "LaplaceTransform.hpp"
#include "LT_interference.hpp"
#include "LT_signal_chi.hpp"
#include "LT_signal_gamma.hpp"
#include "library/GenzMalik/Cube.hpp"
#include "StochasticGeometry/Param.hpp"

class LT_upsilon : public LaplaceTransform {

    LT_interference interference;
    LT_signal_chi signal_chi;
    Param param;
    // LT_signal_gamma signal_gamma;

    public:
    LT_upsilon(const Param& param_, CUBE::Cube<double, 2>& cluster_cube_): interference(param_), signal_chi(param_, cluster_cube_), param(param_){}
    ~LT_upsilon(){}

    std::complex<double> lt_eval(const std::complex<double>& s, const double& threshold, const double& distance) override {
        std::complex<double> interference_eval = interference.lt_eval(s, threshold, distance);
        std::complex<double> signal_eval{0.0, 0.0};
        if (param.RISden != 0.0){
            signal_eval = signal_chi.lt_eval(s, threshold, distance);
        }
        std::complex<double> noise_term = s / PL_d(distance) * param.noise_level / param.Antenna_gain * threshold;
        std::complex<double> eval_result = std::exp(- noise_term - param.BSden * interference_eval - param.RISden * signal_eval); //  

        return eval_result;

    }

};