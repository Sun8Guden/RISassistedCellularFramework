#pragma once

#include <complex>
#include "LaplaceTransform.hpp"
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"
#include "../StochasticGeometryHelper/Param.hpp"
#include "LT_utils.hpp"

class LT_signal_chi : public LaplaceTransform {

    CUBE::Cube<double, 2> cluster_cube;
    Param param;
    
    double zeta_mean {0.821658900384983};
    double zeta_variance {0.324876651418140};
    double s_eff_default {0.0};
    double nc_param {0.0};

    public:
    LT_signal_chi(const Param& param_, CUBE::Cube<double, 2>& cluster_cube_):cluster_cube(cluster_cube_), param(param_){
        s_eff_default = param.Antenna_gain * param.RISnumber / param.UEs_cell * zeta_variance;
        nc_param = (param.RISnumber / param.UEs_cell) * zeta_mean * zeta_mean / zeta_variance;
    }
    ~LT_signal_chi(){}
    
    std::complex<double> lt_eval(const std::complex<double>& s, const double& threshold, const double& distance) override {
        auto f_sig = [&](const double y, const double theta){
            double PL_coef = PL_r1(y) * PL_r2(distance, y, theta) / PL_d(distance);
            double func_sig_eval = PL_coef * s_eff_default;
            std::complex<double> lt_chi = exp(s * nc_param * func_sig_eval / (1.0 - s * 2.0 * func_sig_eval)) \
                /sqrt(1.0 - s * 2.0 * func_sig_eval);
            std::complex<double> ret_val = y * (1.0 - lt_chi);
                    // std::cout << "Test lt_chi: "<< lt_chi << ", and " << y << std::endl;
            return ret_val;
        };
        double esti_error{0.0};
        unsigned int num_eval{0};
        auto cur_cluster_cube = cluster_cube;
        std::complex<double> lt_signal = GM::GM2D<double>::integrate(f_sig, cur_cluster_cube, 1e-6, esti_error, 5000, num_eval);
        return lt_signal;
    }

};