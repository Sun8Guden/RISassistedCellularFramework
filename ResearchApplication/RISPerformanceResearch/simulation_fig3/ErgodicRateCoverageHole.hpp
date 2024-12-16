#pragma once
#include <complex>
#include <iostream>
#include <cmath>
#include <array>
#include "../library/cuhre/Integration.hpp"
#include "../library/cuhre/Region.hpp"
#include "./StochasticGeometry/Param.hpp"


template <typename LT>
double ergodic_rate_fix_distance_coverage_hole(LT lt_interference_signal, std::complex<double> s_parameter, double distance_arg, double x_UE, double y_UE)
{
    double threshold_max {1e4}; double cauchy {40};
    std::array<double, 2> center2d{threshold_max/2.0,  0.0};
    std::array<double, 2> width2d{threshold_max,  cauchy};
    Region<2> cur_region2d = Region<2>(center2d, width2d);
    IntegrationStrategy parallel = IntegrationStrategy::parallel;

    double esti_integral_1{ 0.0 }, esti_integral_2{ 0.0 };
    double esti_error_1{ 0.0 }, esti_error_2{ 0.0 };
    int num_f_eval_1{ 0 }, num_f_eval_2{ 0 };

    auto f_cpv = [&](const std::array<double, 2>& input2d) {
        double z = std::exp(input2d.at(1));
        std::complex<double> z_imag =  std::complex<double>(0.0, 1.0) * z;
        std::complex<double> lt_vec_part_one = lt_interference_signal( s_parameter - z_imag, distance_arg, input2d.at(0), x_UE, y_UE );
        std::complex<double> lt_part_one_zero = lt_interference_signal( - z_imag, distance_arg, input2d.at(0), x_UE, y_UE);
        std::complex<double> cpv =  lt_vec_part_one - lt_part_one_zero;  
        double cpv_val = cpv.imag();
        double ret_val = cpv_val / (1.0 + input2d.at(0) ) / M_PI;
        return ret_val;
    };
    esti_integral_1 = Integration<2>::integrate(f_cpv, cur_region2d, esti_error_1, num_f_eval_1, parallel, 1e-6, 500, 50000);

    // std::cout <<"esti_integral_1 " << esti_integral_1 << std::endl; 

    auto f_part2 = [&](double threshold ){
        std::complex<double> lt_vec_part_one = lt_interference_signal(s_parameter, distance_arg, threshold, x_UE, y_UE );
        std::complex<double> lt_part2 = 1.0 + lt_vec_part_one;
        double part_2_val = 0.5 * lt_part2.real();
        double ret_val = part_2_val / ( 1.0 + threshold );
        return ret_val;
    };
    esti_integral_2 =  gauss_kronrod<double, 15>::integrate(f_part2, 0.0, threshold_max, 10, 1e-6, &esti_error_2); 
    // std::cout <<"esti_integral_2 " << esti_integral_2 << std::endl; 


    double esti_ergodic_rate = esti_integral_1 + esti_integral_2;
    double esti_error = esti_error_1 + esti_error_2;
    unsigned int num_f_eval = num_f_eval_1 + num_f_eval_2; 

    // std::cout << "The integral is " <<  esti_ergodic_rate << " Number of function evaluation: " << \
    //     num_f_eval << ", with total estimated relative error " << esti_error / esti_ergodic_rate << std::endl;

    return esti_ergodic_rate; 
}





// template <typename LT>
// double ergodic_rate_coverage_hole(LT lt_interference_signal, std::complex<double> s_parameter, double distance_arg, double R_hole)
// {
//     double threshold_max {1e4}; double cauchy {40};
//     std::array<double, 4> center4d{threshold_max/2.0,  0.0, R_hole/2.0, M_PI};
//     std::array<double, 4> width4d{threshold_max,  cauchy,  R_hole, 2.0 * M_PI};
//     Region<4> cur_region4d = Region<4>(center4d, width4d);

//     std::array<double, 3> center3d{threshold_max/2.0, R_hole/2.0, M_PI};
//     std::array<double, 3> width3d{threshold_max,  R_hole, 2.0 * M_PI};
//     Region<3> cur_region3d = Region<3>(center3d, width3d);

//     IntegrationStrategy parallel = IntegrationStrategy::parallel;

//     double esti_integral_1{ 0.0 }, esti_integral_2{ 0.0 };
//     double esti_error_1{ 0.0 }, esti_error_2{ 0.0 };
//     int num_f_eval_1{ 0 }, num_f_eval_2{ 0 };





//     auto f_cpv = [&](const std::array<double, 4>& input4d) {
//         double z = std::exp(input4d.at(1));
//         std::complex<double> z_imag =  std::complex<double>(0.0, 1.0) * z;
//         std::complex<double> lt_vec_part_one = lt_interference_signal( s_parameter - z_imag, distance_arg, input4d.at(0), input4d.at(2), input4d.at(3));
//         std::complex<double> lt_part_one_zero = lt_interference_signal( - z_imag, distance_arg, input4d.at(0), input4d.at(2), input4d.at(3) );
//         std::complex<double> cpv =  lt_vec_part_one - lt_part_one_zero;  
//         double cpv_val = cpv.imag();
//         double ret_val = cpv_val / (1.0 + input4d.at(0) ) / M_PI;
//         return ret_val; // * input4d.at(2);
//     };
//     esti_integral_1 = Integration<4>::integrate(f_cpv, cur_region4d, esti_error_1, num_f_eval_1, parallel, 1e-6, 500, 500000);

//     std::cout <<"esti_integral_1 " << esti_integral_1 << std::endl; 

//     auto f_part2 = [&](const std::array<double, 3>& input3d ){

//         std::complex<double> lt_vec_part_one = lt_interference_signal(s_parameter, distance_arg, input3d.at(0), input3d.at(1), input3d.at(2) );
//         std::complex<double> lt_part2 = 1.0 + lt_vec_part_one;
//         double part_2_val = 0.5 * lt_part2.real();
//         double ret_val = part_2_val / ( 1.0 + input3d.at(0) );
//         return ret_val; // * input3d.at(1);
//     };
//     esti_integral_2 =   Integration<3>::integrate(f_part2, cur_region3d, esti_error_2, num_f_eval_2, parallel, 1e-6, 500, 50000);
//     std::cout <<"esti_integral_2 " << esti_integral_2 << std::endl; 


//     double esti_ergodic_rate = esti_integral_1 + esti_integral_2;
//     double esti_error = esti_error_1 + esti_error_2;
//     unsigned int num_f_eval = num_f_eval_1 + num_f_eval_2; 

//     // std::cout << "The integral is " <<  esti_ergodic_rate << " Number of function evaluation: " << \
//     //     num_f_eval << ", with total estimated relative error " << esti_error / esti_ergodic_rate << std::endl;

//     return esti_ergodic_rate / (M_PI * R_hole *R_hole); 
// }