#pragma once

#include <cmath>
#include <complex>
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../LaplaceTransformInstance/Interference.hpp"
#include "../StochasticGeometry/Param.hpp"



template <typename S_type>
S_type LT_interference_mean_BS(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP){

    S_type exp_part = LT_interference_BSs(s, distance, threshold, param, BS_PP);
    S_type no_exp_part = LT_interference_BSs_D1_no_exp(s, distance, threshold, param, BS_PP);
    return no_exp_part * exp_part;

}


template <typename S_type>
S_type LT_interference_mean_scatter(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP, PoissonPP& ScatteredInterferencePP){
    
    S_type RIS_scatter_exp_part = LT_interference_RIS_scatter(s, distance, threshold, param, BS_PP, ScatteredInterferencePP);
    S_type RIS_scatter_no_exp_part = LT_interference_RIS_scatter_D1_no_exp(s, distance, threshold, param, BS_PP, ScatteredInterferencePP);

    return RIS_scatter_exp_part * RIS_scatter_no_exp_part;
}

template <typename S_type>
S_type LT_interference_mean_beamform(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP, PoissonPP& BeamformedInterferencePP){
    
    S_type RIS_beam_exp_part = LT_interference_RIS_beam(s, distance, threshold, param, BS_PP, BeamformedInterferencePP);
    S_type RIS_beam_no_exp_part = LT_interference_RIS_beam_D1_no_exp(s, distance, threshold, param, BS_PP, BeamformedInterferencePP);

    return RIS_beam_exp_part * RIS_beam_no_exp_part;
}

template <typename S_type>
S_type LT_interference_mean_RIS(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    S_type RIS_exp_part = LT_interference_RIS(s, distance, threshold, param, BS_PP, BeamformedInterferencePP, ScatteredInterferencePP);
    S_type RIS_no_exp_part = LT_interference_RIS_D1_no_exp(s, distance, threshold, param, BS_PP, BeamformedInterferencePP, ScatteredInterferencePP);

    return RIS_exp_part * RIS_no_exp_part;
}

template <typename S_type>
S_type LT_interference_mean_BS_RIS(S_type s, double distance, double threshold, const Param& param, PoissonPP& BS_PP,\
         PoissonPP& BeamformedInterferencePP, PoissonPP& ScatteredInterferencePP){

    S_type BS_RIS_exp_part = LT_interference(s, distance, threshold, param, BS_PP, BeamformedInterferencePP, ScatteredInterferencePP);
    S_type BS_RIS_no_exp_part = LT_interference_D1_no_exp(s, distance, threshold, param, BS_PP, BeamformedInterferencePP, ScatteredInterferencePP);

    return BS_RIS_exp_part * BS_RIS_no_exp_part;
}

