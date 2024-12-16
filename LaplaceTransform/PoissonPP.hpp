#pragma once
#include <iostream>
#include <cassert>
#include <complex>
#include <utility>
#include <limits>
#include <functional>
#include <cmath>
#include "LaplaceTransform.hpp"
#include "../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#include "../Integration/Cuhre/Integration.hpp"
#include "../Integration/Cuhre/Region.hpp"
#include "../Integration/GenzMalik/Cube.hpp"
#include "../Integration/GenzMalik/GM2D.hpp"

using namespace boost::math::quadrature;

// We introduce two classes for 2D Poisson point processes: one for centered scenario 
// where we can reduce one dimension whereas the other for non-centered scenario. 
// They are mainly interest for stochastic geometry research.

// We also provide the requirement for an interface for general implementation.



class PoissonPP2D_Reduced: public LaplaceTransform {
    // 
    double m_density_pp;
    std::function<std::complex<double>(std::complex<double>, double)> mf_laplace_function_1D;
    double m_region_left;
    double m_region_right;

    public:
    PoissonPP2D_Reduced(const std::function<std::complex<double>(std::complex<double>, double)>& \
        laplace_function_1D, double density_pp, double x_left, double x_right = std::numeric_limits<double>::infinity()\
        ):mf_laplace_function_1D(laplace_function_1D), m_density_pp(density_pp), m_region_left(x_left), m_region_right(x_right){};
    
    void set_functor(const std::function<std::complex<double>(std::complex<double>, double)>& laplace_function_1D){
        mf_laplace_function_1D = laplace_function_1D;
    };

    void set_region(double x_left, double x_right = std::numeric_limits<double>::infinity() ){
        m_region_left = x_left;
        m_region_right = x_right;
    };

    std::complex<double> eval(const std::complex<double>& s) const override {
        double estimated_error;
        auto f_integrand = [&](const double distance){
            return distance * (1.0 - mf_laplace_function_1D(s, distance));
        };
        std::complex<double> integral_value = gauss_kronrod<double, 15>::integrate(f_integrand, \
        m_region_left,  m_region_right, 10, 1e-6, &estimated_error);
        return std::exp( - 2.0 * M_PI * integral_value * m_density_pp);
    };

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    std::complex<double> eval_exponent(const std::complex<double>& s) const override {
        double estimated_error;
        auto f_integrand = [&](const double distance){
            return distance * (1.0 - mf_laplace_function_1D(s, distance));
        };
        std::complex<double> integral_value = gauss_kronrod<double, 15>::integrate(f_integrand, \
        m_region_left,  m_region_right, 10, 1e-6, &estimated_error);
        return  - 2.0 * M_PI * integral_value;
    };

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};


class PoissonPP2D : public LaplaceTransform {
        std::function<std::complex<double>(std::complex<double>, double, double)> mf_laplace_transform_2D;
        CUBE::Cube<double, 2> m_cluster_cube; 
        double m_density_pp; 
    public:
        PoissonPP2D(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D, \
            double density_pp, const CUBE::Cube<double, 2>& cluster_cube): mf_laplace_transform_2D(laplace_function_2D),\
            m_density_pp(density_pp), m_cluster_cube(cluster_cube){}

        void set_functor(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D){
            mf_laplace_transform_2D = laplace_function_2D;
        }

        void set_region(const CUBE::Cube<double, 2>& cluster_cube){
            m_cluster_cube = cluster_cube;
        }

    std::complex<double> eval(const std::complex<double>& s) const override {

        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        auto f_integrand = [&](const double y, const double theta){
            return y * (1.0 - mf_laplace_transform_2D(s, y, theta));
        };

        std::complex<double> integral_value = GM::GM2D<double>::integrate(f_integrand, \
            m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return std::exp( - integral_value * m_density_pp );
    }

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    std::complex<double> eval_exponent(const std::complex<double>& s) const override {

        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        auto f_integrand = [&](const double y, const double theta){
            return y * (1.0 - mf_laplace_transform_2D(s, y, theta));
        };

        std::complex<double> integral_value = GM::GM2D<double>::integrate(f_integrand, \
            m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return  - integral_value;
    }

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};

