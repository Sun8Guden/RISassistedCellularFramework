#pragma once
#include <cmath>
#include <complex>
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"
#include "../LaplaceTransformDistribution/ChiSquare.hpp"
#include "../StochasticGeometryHelper/Param.hpp"
#include "../LaplaceTransformFunctions/LT_utils.hpp"



template <typename S_type>
S_type LaplaceReflectedSignalsWedgePPP(S_type s, double distance, const CUBE::Cube<double, 2>& cluster_cube, const Param& param){

    auto f_reflected_signal = [&](double y, double theta){

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_r2(distance, y, theta) / PL_d(distance);
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * (1.0+1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * (1.0-1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        
        S_type return_value =  y * (1.0 - ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) * ChiSquare::eval_central(s, central_variance));
        return return_value;
    };

    double estimated_error{0.0};
    unsigned int number_evaluation{0};
    auto integration_value = GM::GM2D<double>::integrate(f_reflected_signal, cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);

    S_type return_value = std::exp( - param.RISden * integration_value);
    return return_value;
}


template <typename S_type>
S_type LaplaceReflectedSignalsWedgeBPP(S_type s, double distance, const CUBE::Cube<double, 2>& cluster_cube, const Param& param, double number_points){

    auto f_reflected_signal = [&](double y, double theta){

        double PL_coef = param.Antenna_gain * PL_r1(y) * PL_r2(distance, y, theta) / PL_d(distance);
        double non_central_mean = PL_coef * (param.RISnumber / param.UEs_cell) * (param.RISnumber / param.UEs_cell) * param.zeta_mean * param.zeta_mean;
        double non_central_variance =  PL_coef * (1.0+1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        double central_variance =  PL_coef * (1.0-1.0/param.UEs_cell)/2.0 * param.RISnumber * param.zeta_variance;
        
        S_type return_value =  y * ChiSquare::eval_non_central(s, non_central_mean, non_central_variance) * ChiSquare::eval_central(s, central_variance);
        return return_value;
    };

    double estimated_error{0.0};
    unsigned int number_evaluation{0};
    auto integration_value = GM::GM2D<double>::integrate(f_reflected_signal, cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);

    double radius_max = cluster_cube.cube_center.at(0) + 0.5 * cluster_cube.cube_size.at(0);
    double radius_min = cluster_cube.cube_center.at(0) - 0.5 * cluster_cube.cube_size.at(0);
    double angle_range = cluster_cube.cube_size.at(1);
    double BPP_probability =  ( angle_range * 0.5 * (radius_max * radius_max - radius_min * radius_min) );


    S_type return_value = std::pow(( integration_value / BPP_probability ), number_points);
    return return_value;
}
