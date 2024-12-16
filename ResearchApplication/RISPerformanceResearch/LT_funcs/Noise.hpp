#pragma once
#include <complex>
#include <cmath>
#include "../StochasticGeometry/Param.hpp"
#include "LT_utils.hpp"

template <typename S_type>
S_type LT_noise(S_type s, double noise_power, double distance, double threshold, const Param& param){

    double direct_interference_variance = threshold * noise_power / (param.Antenna_gain * PL_d(distance));
    S_type return_value = std::exp( - s * direct_interference_variance);

    return return_value;
}