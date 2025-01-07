#pragma once
#include <complex>
#include <cmath>
#include "../LaplaceTransformDistribution/Exponential.hpp"
#include "../LaplaceTransformDistribution/ChiSquare.hpp"
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../StochasticGeometryHelper/Param.hpp"
#include "LT_utils.hpp"
#include "ReflectedSignal.hpp"


template <typename S_type>
S_type LaplaceReflectedInterference(S_type s, double cur_BS_distance, double distance, double threshold, const Param& param, PoissonPP& PPP){

    auto f_reflected_signal = [&](double y, double theta){
        // For signal, s is negative. (probability I will add an assertion here.)

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_dN(cur_BS_distance, y, theta) / PL_d(distance) * threshold;
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * ( 1.0 + 1.0 / param.UEs_cell ) / 2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * ( 1.0 - 1.0 / param.UEs_cell ) / 2.0 * param.RISnumber * param.zeta_variance;

        S_type return_value = y * (1.0 - ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) * ChiSquare::eval_central(s, central_variance));
        return return_value;
    };

    S_type return_value = PPP.eval_2D(f_reflected_signal);
    return return_value;
}

template <typename S_type>
S_type LaplaceReflectedInterference_D1_no_exp(S_type s, double cur_BS_distance, double distance, double threshold, const Param& param, PoissonPP& PPP){

    auto f_reflected_signal = [&](double y, double theta){
        // For signal, s is negative. (probability I will add an assertion here.)

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_dN(cur_BS_distance, y, theta) / PL_d(distance) * threshold;
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * ( 1.0 + 1.0 / param.UEs_cell ) / 2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * ( 1.0 - 1.0 / param.UEs_cell ) / 2.0 * param.RISnumber * param.zeta_variance;

        S_type return_value = y * ( - ChiSquare::eval_non_central_D1(s, non_central_mean, non_central_variance) * ChiSquare::eval_central(s, central_variance)\
         - ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) * ChiSquare::eval_central_D1(s, central_variance));
        return return_value;
    };

    S_type return_value = PPP.eval_2D_no_exp(f_reflected_signal);
    return return_value;
}

template <typename S_type>
S_type LaplaceScatteredInterference(S_type s,  double cur_BS_distance, double distance, double threshold, const Param& param, PoissonPP& ScatteredInterferencePP){

    auto f_reflected_signal = [&](double y, double theta){

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_dN(cur_BS_distance, y, theta) / PL_d(distance) * threshold;
        double central_variance = PL_coef * ( 1.0 - 1.0 / param.UEs_cell ) * param.RISnumber * param.zeta_variance;
        return y * ( 1.0 - Exponential::eval(s, central_variance));

    };

    S_type return_value = ScatteredInterferencePP.eval_2D(f_reflected_signal);
    return return_value;
}

template <typename S_type>
S_type LaplaceScatteredInterference_D1_no_exp(S_type s,  double cur_BS_distance, double distance, double threshold, const Param& param, PoissonPP& ScatteredInterferencePP){

    auto f_reflected_signal = [&](double y, double theta){

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_dN(cur_BS_distance, y, theta) / PL_d(distance) * threshold;
        double central_variance = PL_coef * ( 1.0 - 1.0 / param.UEs_cell ) * param.RISnumber * param.zeta_variance;
        return y * ( - Exponential::eval_D1(s, central_variance));

    };

    S_type return_value = ScatteredInterferencePP.eval_2D_no_exp(f_reflected_signal);
    return return_value;
}

template <typename S_type>
S_type LT_interference_BSs(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP){

    auto f_interference = [&](double x){

        double direct_interference_variance = threshold * PL_d(x) / PL_d(distance);
        S_type direct_interference_value = Exponential::eval(s, direct_interference_variance);
        S_type return_value = x * ( 1.0 - direct_interference_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D(f_interference, distance, 3000.0);

    return return_value;
}

template <typename S_type>
S_type LT_interference_BSs_D1_no_exp(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP){

    auto f_interference = [&](double x){
        double direct_interference_variance = threshold * PL_d(x) / PL_d(distance);
        S_type return_value = Exponential::eval_D1(s, direct_interference_variance);
        return_value = x * (- return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D_no_exp(f_interference, distance, 3000.0);

    return return_value;
}


template <typename S_type>
S_type LT_interference_RIS(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){
        S_type scattered_interference = LaplaceScatteredInterference(s, x, distance, threshold, param, ScatteredInterferencePP);
        S_type reflected_interference = LaplaceReflectedInterference(s, x, distance, threshold, param, BeamformedInterferencePP);
        S_type return_value =  reflected_interference * scattered_interference;
        return_value = x * ( 1.0 - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D(f_interference, distance, 3000.0);
    return return_value;
}


template <typename S_type>
S_type LT_interference_RIS_beam(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
        PoissonPP& BeamformedInterferencePP){

    auto f_interference = [&](double x){
        S_type reflected_interference = LaplaceReflectedInterference(s, x, distance, threshold, param, BeamformedInterferencePP);
        S_type return_value =  reflected_interference;
        return_value = x * ( 1.0 - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D(f_interference, distance, 3000.0);
    return return_value;
}


template <typename S_type>
S_type LT_interference_RIS_scatter(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
        PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){
        S_type scattered_interference = LaplaceScatteredInterference(s, x, distance, threshold, param, ScatteredInterferencePP);
        S_type return_value = scattered_interference;
        return_value = x * ( 1.0 - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D(f_interference, distance, 3000.0);
    return return_value;
}

template <typename S_type>
S_type LT_interference_RIS_D1_no_exp(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){
        S_type scattered_interference = LaplaceScatteredInterference(s, x, distance, threshold, param, ScatteredInterferencePP);
        S_type scattered_interference_D1 = LaplaceScatteredInterference_D1_no_exp(s, x, distance, threshold, param, ScatteredInterferencePP);

        S_type reflected_interference = LaplaceReflectedInterference(s, x, distance, threshold, param, BeamformedInterferencePP);
        S_type reflected_interference_D1 = LaplaceReflectedInterference_D1_no_exp(s, x, distance, threshold, param, BeamformedInterferencePP);

        S_type return_value = reflected_interference * scattered_interference *(scattered_interference_D1 + reflected_interference_D1);
        return_value = x * ( - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D_no_exp(f_interference, distance, 3000.0);

    return return_value;
}

template <typename S_type>
S_type LT_interference_RIS_beam_D1_no_exp(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP){

    auto f_interference = [&](double x){

        S_type return_value = LaplaceReflectedInterference_D1_no_exp(s, x, distance, threshold, param, BeamformedInterferencePP);
        return_value = x * ( - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D_no_exp(f_interference, distance, 3000.0);

    return return_value;
}

template <typename S_type>
S_type LT_interference_RIS_scatter_D1_no_exp(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){
        S_type return_value = LaplaceScatteredInterference_D1_no_exp(s, x, distance, threshold, param, ScatteredInterferencePP);
        return_value = x * ( - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D_no_exp(f_interference, distance, 3000.0);

    return return_value;
}


template <typename S_type>
S_type LT_interference(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){

        double direct_interference_variance = threshold * PL_d(x) / PL_d(distance);
        S_type direct_interference_value = Exponential::eval(s, direct_interference_variance);
        S_type scattered_interference = LaplaceScatteredInterference(s, x, distance, threshold, param, ScatteredInterferencePP);
        S_type reflected_interference = LaplaceReflectedInterference(s, x, distance, threshold, param, BeamformedInterferencePP);
        S_type return_value = direct_interference_value * reflected_interference * scattered_interference;
        return_value = x * ( 1.0 - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D(f_interference, distance, 3000.0);

    return return_value;
}

template <typename S_type>
S_type LT_interference_D1_no_exp(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    auto f_interference = [&](double x){

        double direct_interference_variance = threshold * PL_d(x) / PL_d(distance);
        S_type direct_interference_value = Exponential::eval(s, direct_interference_variance);
        S_type direct_interference_D1_value = Exponential::eval_D1(s, direct_interference_variance);

        S_type scattered_interference = LaplaceScatteredInterference(s, x, distance, threshold, param, ScatteredInterferencePP);
        S_type scattered_interference_D1 = LaplaceScatteredInterference_D1_no_exp(s, x, distance, threshold, param, ScatteredInterferencePP);

        S_type reflected_interference = LaplaceReflectedInterference(s, x, distance, threshold, param, BeamformedInterferencePP);
        S_type reflected_interference_D1 = LaplaceReflectedInterference_D1_no_exp(s, x, distance, threshold, param, BeamformedInterferencePP);

        S_type return_value = reflected_interference * scattered_interference *( direct_interference_D1_value  +\
                            direct_interference_value * (scattered_interference_D1 + reflected_interference_D1));
        return_value = x * ( - return_value );
        return return_value;
    };
    
    S_type return_value = BS_PP.eval_1D_no_exp(f_interference, distance, 3000.0);

    return return_value;
}

