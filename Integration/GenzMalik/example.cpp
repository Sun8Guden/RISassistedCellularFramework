#include "GenzMalik/GM2D.hpp"
#include "GenzMalik/GM3D.hpp"
#include "GenzMalik/Cube.hpp"
#include <iostream>
#include <array>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iomanip>

// This is a trivial example simply to showcase how to use 

int main(){

    double estimated_error {0.0};

    CUBE::Cube<double, 3> cube3 = CUBE::make_cube_3D<double>(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

    auto foo_GM = [&](const double& x, const double& y, const double& z){
        double arg_value = - (1.0+x)/(1.0-x)  - 3.0 * (1.0+y)/(1.0-y) - (1.0+z)/(1.0-z)*(1.0+z)/(1.0-z);
        double exp_value = std::exp(arg_value);
        double ret_value = exp_value * 8.0 / (1.0-x)/(1.0-x)/(1.0-y)/(1.0-y)/(1.0-z)/(1.0-z);
        return ret_value; 
    };
    unsigned int num_func_eval2{0};
    auto integral_2 =  GM::GM3D<double>::integrate(foo_GM, cube3, 1e-6, estimated_error, 30,  num_func_eval2); 
    // std::cout << std::setprecision(15)<< "The numerical integral using Genz Malik method is " << integral_2 << ", with an estimated error " << estimated_error << " after " << num_func_eval2 << " function calls\n";
    return 0;
}