#pragma once
#include <cmath>
#include <complex>
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../LaplaceTransformPointProcess/BinomialPP.hpp"
#include "../LaplaceTransformDistribution/ChiSquare.hpp"
#include "../StochasticGeometryHelper/Param.hpp"
#include "LT_utils.hpp"


template <typename S_type>
S_type LaplaceReflectedSignals(S_type s, double distance, const Param& param, PoissonPP& PPP){

    auto f_reflected_signal = [&](double y, double theta){
        // For signal, s is negative. (probability I will add an assertion here.)

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_r2(distance, y, theta) / PL_d(distance);
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * (1.0+1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * (1.0-1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        
        S_type return_value = y * (1.0 - ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) *ChiSquare::eval_central(s, central_variance));
        return return_value;
    };
 
    S_type return_value = PPP.eval_2D(f_reflected_signal);
    return return_value;
}


template <typename S_type>
S_type LaplaceReflectedSignals(S_type s, double distance, const Param& param, BinomialPP& BPP){

    auto f_reflected_signal = [&](double y, double theta){

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_r2(distance, y, theta) / PL_d(distance);
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * (1.0+1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * (1.0-1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        
        S_type return_value = y * ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) *ChiSquare::eval_central(s, central_variance);
        // std::cout << "return_value " << return_value << std::endl;
        return return_value;
    };
 
    S_type return_value = BPP.eval_2D(f_reflected_signal);
    return return_value;
}