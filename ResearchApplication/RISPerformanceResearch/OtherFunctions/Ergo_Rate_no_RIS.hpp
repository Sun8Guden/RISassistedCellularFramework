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
#include "boost/math/quadrature/gauss_kronrod.hpp"
#include "Param.hpp"
#include "Cov_Prob.hpp"

class Ergo_Rate_no_RIS {
    // The third argument will be exponential of the value, which are infinitely small or large. 

    // CUBE::Cube<double, 2> cube2 = CUBE::make_cube_2D<double>(double(0.0), double(5000.0), double(30.0), double(2000.0) );

    // exp_sinh<double> singular;
    Region<2> cur_region2d;
    IntegrationStrategy parallel = IntegrationStrategy::parallel;
    Param param;

public:
    Ergo_Rate_no_RIS( const Param& param_ ) : param( param_ ) {
        std::array<double, 4> region_info2d{0.0, 5000.0, 30.0, 2000.0};
        std::array<double, 2> center2d{(region_info2d.at(0)+region_info2d.at(1))/2.0, (region_info2d.at(2)+region_info2d.at(3))/2.0};
        std::array<double, 2> width2d{region_info2d.at(1)-region_info2d.at(0), region_info2d.at(3)-region_info2d.at(2)};
        cur_region2d = Region<2>(center2d, width2d);
    }
    ~Ergo_Rate_no_RIS(){}
    
    double calc(LaplaceTransform& lt_trans, const double& s_BF) {
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

};