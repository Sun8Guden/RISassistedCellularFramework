#pragma once
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <complex>
// #include "library/GenzMalik/GM3D.hpp"
// #include "library/GenzMalik/GM2D.hpp"
#include "library/cuhre/Integration.hpp"
#include "library/cuhre/Region.hpp"
#include "LaplacetransformArchiv/LaplaceTransform.hpp"
#include "Param.hpp"
#include "Cov_Prob.hpp"

class Ergo_Rate_BF {
    // The third argument will be exponential of the value, which are infinitely small or large. 

    Region<2> cur_region2d;
    Region<3> cur_region3d;
    IntegrationStrategy parallel = IntegrationStrategy::parallel;

    // CUBE::Cube<double, 3> cube3 = CUBE::make_cube_3D<double>(double(0.0), double(5000.0), double(30.0), double(2000.0), double(-20.0), double(20.0) );    
    // CUBE::Cube<double, 2> cube2 = CUBE::make_cube_2D<double>(double(0.0), double(5000.0), double(30.0), double(2000.0) );
    // exp_sinh<double> singular;

    std::complex<double> imag_unit;
    Param param;

public:
    Ergo_Rate_BF( const Param& param_ ) : imag_unit(0.0, 1.0), param( param_ ) {
        std::array<double, 4> region_info2d{0.0, 5000.0, 30.0, 2000.0};
        std::array<double, 2> center2d{(region_info2d.at(0)+region_info2d.at(1))/2.0, (region_info2d.at(2)+region_info2d.at(3))/2.0};
        std::array<double, 2> width2d{region_info2d.at(1)-region_info2d.at(0), region_info2d.at(3)-region_info2d.at(2)};
        cur_region2d = Region<2>(center2d, width2d);

        std::array<double, 6> region_info3d{0.0, 5000.0, 30.0, 2000.0, -20.0, 20.0};
        std::array<double, 3> center3d{(region_info3d.at(0)+region_info3d.at(1))/2.0, (region_info3d.at(2)+region_info3d.at(3))/2.0, (region_info3d.at(4)+region_info3d.at(5))/2.0 };
        std::array<double, 3> width3d{region_info3d.at(1)-region_info3d.at(0), region_info3d.at(3)-region_info3d.at(2),  region_info3d.at(5)-region_info3d.at(4)};
        cur_region3d = Region<3>(center3d, width3d);
    }
    ~Ergo_Rate_BF(){}

    double calc_no_RIS(LaplaceTransform& lt_trans, const double& s_BF) {
        double  esti_error_2{ 0.0 };
        double  esti_integral_2{ 0.0 };
        int  num_f_eval_2{ 0 };

        auto f_part2 = [&](const std::array<double, 2>& input2d){
            std::complex<double> lt_vec_part_one = lt_trans.lt_eval(s_BF, input2d.at(0), input2d.at(1));
            double part_2_val = lt_vec_part_one.real();
            double ret_val = part_2_val / ( 1.0 + input2d.at(0) ) * std::exp(- M_PI * param.BSden * input2d.at(1) * input2d.at(1) ) * 2.0 * M_PI * input2d.at(1) * param.BSden;
            return ret_val;
        };
        esti_integral_2 =  Integration<2>::integrate(f_part2, cur_region2d, esti_error_2, num_f_eval_2, parallel, 1e-6, 500, 100000); 

        std::cout << "The integral is " <<  esti_integral_2 << " Number of function evaluation: " << \
            num_f_eval_2 << ", with estimated relative error " << esti_error_2 << std::endl;

        return esti_integral_2; 
    }
    
    double calc(LaplaceTransform& lt_trans, const double& s_BF) {
        double esti_integral_1{ 0.0 }, esti_integral_2{ 0.0 };
        double esti_error_1{ 0.0 }, esti_error_2{ 0.0 };
        int num_f_eval_1{ 0 }, num_f_eval_2{ 0 };

        auto f_cpv = [&](const std::array<double, 3>& input3d) {
            double z = exp(input3d.at(2));
            std::complex<double> z_imag =  imag_unit * z;

            std::complex<double> lt_vec_part_one = lt_trans.lt_eval( s_BF - z_imag, input3d.at(0), input3d.at(1) );
            std::complex<double> lt_part_one_zero = lt_trans.lt_eval( - z_imag, input3d.at(0), input3d.at(1) );
            std::complex<double> cpv =  lt_vec_part_one - lt_part_one_zero;  
            double cpv_val = cpv.imag();
            double ret_val = cpv_val / (1.0 + input3d.at(0)) * exp(- M_PI * param.BSden * input3d.at(1) * input3d.at(1) ) * 2.0 * input3d.at(1) * param.BSden;
            return ret_val;
        };
        esti_integral_1 = Integration<3>::integrate(f_cpv, cur_region3d, esti_error_1, num_f_eval_1, parallel, 1e-6, 500, 400000);


        auto f_part2 = [&](const std::array<double, 2>& input2d){
            std::complex<double> lt_vec_part_one = lt_trans.lt_eval(s_BF, input2d.at(0), input2d.at(1) );
            std::complex<double> lt_part2 = 1.0 + lt_vec_part_one;
            double part_2_val = 0.5 * lt_part2.real();
            double ret_val = part_2_val / ( 1.0 + input2d.at(0) ) * exp(- M_PI * param.BSden * input2d.at(1) * input2d.at(1) ) * 2.0 * M_PI * input2d.at(1) * param.BSden;
            return ret_val;
        };
        esti_integral_2 =  Integration<2>::integrate(f_part2, cur_region2d, esti_error_2, num_f_eval_2, parallel, 1e-6, 500, 20000); 
        double esti_ergodic_rate = esti_integral_1 + esti_integral_2;
        double esti_error = esti_error_1 + esti_error_2;
        unsigned int num_f_eval = num_f_eval_1 + num_f_eval_2; 

        // std::cout << "The integral is " <<  esti_ergodic_rate << " Number of function evaluation: " << \
        //     num_f_eval << ", with total estimated relative error " << esti_error / esti_ergodic_rate << std::endl;

        return esti_ergodic_rate; 
    }

};