#pragma once
#include <iostream>
#include <cassert>
#include <complex>
#include <utility>
#include <cmath>
#include <functional>
#include "LaplaceTransform.hpp"
#include "../Integration/boost/math/quadrature/gauss_kronrod.hpp"
#include "../Integration/Cuhre/Integration.hpp"
#include "../Integration/Cuhre/Region.hpp"
#include "../Integration/GenzMalik/Cube.hpp"
#include "../Integration/GenzMalik/GM2D.hpp"

// Since it is only possible to model a finite number of nodes over a finite area, we only consider the 2D case

class BinomialPP : public LaplaceTransform {

    std::function<std::complex<double>(std::complex<double>, double, double)> mf_laplace_transform_2D;
    CUBE::Cube<double, 2> m_cluster_cube; 
    double m_number_pp; 
    double m_probability;

public:

    BinomialPP(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D, \
            double number_pp, const CUBE::Cube<double, 2>& cluster_cube): mf_laplace_transform_2D(laplace_function_2D),\
            m_number_pp(number_pp), m_cluster_cube(cluster_cube){
        double radius_max = cluster_cube.cube_center.at(0) + 0.5 * cluster_cube.cube_size.at(0);
        double radius_min = cluster_cube.cube_center.at(0) - 0.5 * cluster_cube.cube_size.at(0);
        m_probability = 1.0 / ( M_PI * (radius_max * radius_max - radius_min * radius_min) );
    }
    ~BinomialPP(){}

    void set_functor(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D){
        mf_laplace_transform_2D = laplace_function_2D;
    }

    void set_region(const CUBE::Cube<double, 2>& cluster_cube){
        m_cluster_cube = cluster_cube;
    }

    void set_number_pp(double number_pp){
        m_number_pp=number_pp;
    }

    std::complex<double> eval(const std::complex<double>& s) const override {

        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        std::complex<double> integral_value = GM::GM2D<double>::integrate(std::bind(mf_laplace_transform_2D, s, std::placeholders::_1, std::placeholders::_2), \
            m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return std::pow( m_probability * integral_value,  double( m_number_pp ) );;
    }

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};