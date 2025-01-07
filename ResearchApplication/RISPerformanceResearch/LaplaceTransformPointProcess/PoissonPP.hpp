#pragma once
#include <iostream>
#include <cassert>
#include <complex>
#include <utility>
#include <cmath>
#include "../../../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#include "../../../Integration/Cuhre/Integration.hpp"
#include "../../../Integration/Cuhre/Region.hpp"
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"
// #include "library/GenzMalik/GM3D.hpp"

using namespace boost::math::quadrature;

class PoissonPP {
    private:
        double density;
        CUBE::Cube<double, 2> cluster_cube;
    public:
    PoissonPP(double density_, const CUBE::Cube<double, 2>& cluster_cube_): density(density_), cluster_cube(cluster_cube_) {}
    PoissonPP(double density_):density(density_){}
    ~PoissonPP(){}


    template<class F>
    auto eval_1D(F f, double x_left, double x_right)->decltype(std::declval<F>()(std::declval<double>())) const
    {
        typedef decltype(std::declval<F>()(std::declval<double>())) return_type;
        double estimated_error;
        return_type integral_value = gauss_kronrod<double, 15>::integrate(f, x_left, x_right, 10, 1e-6, &estimated_error);
        integral_value = 2.0 * M_PI * integral_value;
        

        return std::exp( - integral_value * density);
    }

    template<class F>
    auto eval_1D_no_exp(F f, double x_left, double x_right)->decltype(std::declval<F>()(std::declval<double>())) const
    {
        typedef decltype(std::declval<F>()(std::declval<double>())) return_type;
        double estimated_error;
        return_type integral_value = gauss_kronrod<double, 15>::integrate(f, x_left, x_right, 10, 1e-6, &estimated_error);
        integral_value = 2.0 * M_PI * integral_value;
        

        return - integral_value * density;
    }

    template<class F>
    auto eval_2D(F f)->decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) const
    {
        typedef decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) return_type;
        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        return_type integral_value = GM::GM2D<double>::integrate(f, cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return std::exp( - integral_value * density );
    }

    template<class F>
    auto eval_2D_no_exp(F f)->decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) const
    {
        typedef decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) return_type;
        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        return_type integral_value = GM::GM2D<double>::integrate(f, cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return - integral_value * density ;
    }

    
};