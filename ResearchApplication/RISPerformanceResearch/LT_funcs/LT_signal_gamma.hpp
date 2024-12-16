#pragma once

#include <complex>
#include "LaplaceTransform.hpp"
#include "library/GenzMalik/Cube.hpp"
#include "library/GenzMalik/GM2D.hpp"
#include "StochasticGeometryArchiv/Param.hpp"


class LT_signal_gamma : public LaplaceTransform {
    CUBE::Cube<double, 2> cluster_cube;
    Param param;

    public:
    LT_signal_gamma(const Param& param_, CUBE::Cube<double, 2>& cluster_cube_):cluster_cube(cluster_cube), param(param_){
       
    }
    ~LT_signal_gamma(){}
    
    std::complex<double> lt_eval(const std::complex<double>& s, const double& threshold, const double& distance) override {
        auto f_sig = [&](const double y, const double theta){
            std::complex<double> ret_val;
            return ret_val;
        };
        double esti_error{0.0};
        unsigned int num_eval{0};
        auto cur_cluster_cube = cluster_cube;
        std::complex<double> signal_estimated = GM::GM2D<double>::integrate(f_sig, cur_cluster_cube, 1e-6, esti_error, 5000, num_eval);
        std::complex<double> lt_signal = -param.RISden * signal_estimated;
        return lt_signal;
    }

};