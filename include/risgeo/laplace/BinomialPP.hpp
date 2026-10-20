#pragma once
#include <complex>
#include <cmath>
#include <functional>
#include "LaplaceTransform.hpp"
#include "../integration/genz_malik/Cube.hpp"
#include "../integration/genz_malik/GM2D.hpp"

namespace risgeo {

// Since it is only possible to model a finite number of nodes over a finite area, we only consider the 2D case.
// The N points are i.i.d. and uniform over the ring (y, theta) described by the cube, so that
// L(s) = ( p * \int\int y L(s; y, theta) dy dtheta )^N with p = 1 / (pi * (R_max^2 - R_min^2)).
// As for PoissonPP2D, the functor gives the Laplace transform of a single point at (y, theta).

class BinomialPP : public LaplaceTransform {

    std::function<std::complex<double>(std::complex<double>, double, double)> mf_laplace_transform_2D;
    CUBE::Cube<double, 2> m_cluster_cube;
    double m_number_pp;
    double m_probability;

    std::complex<double> expected_single_point(const std::complex<double>& s) const {
        double estimated_error{0.0};
        unsigned int number_evaluation{0};
        auto f_integrand = [&](const double y, const double theta){
            return y * mf_laplace_transform_2D(s, y, theta);
        };
        std::complex<double> integral_value = GM::GM2D<double>::integrate(f_integrand, \
            m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        return m_probability * integral_value;
    }

    void update_probability(){
        double radius_max = m_cluster_cube.cube_center.at(0) + 0.5 * m_cluster_cube.cube_size.at(0);
        double radius_min = m_cluster_cube.cube_center.at(0) - 0.5 * m_cluster_cube.cube_size.at(0);
        m_probability = 1.0 / ( M_PI * (radius_max * radius_max - radius_min * radius_min) );
    }

public:

    BinomialPP(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D, \
            double number_pp, const CUBE::Cube<double, 2>& cluster_cube): mf_laplace_transform_2D(laplace_function_2D),\
            m_cluster_cube(cluster_cube), m_number_pp(number_pp){
        update_probability();
    }
    ~BinomialPP(){}

    void set_functor(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D){
        mf_laplace_transform_2D = laplace_function_2D;
    }

    void set_region(const CUBE::Cube<double, 2>& cluster_cube){
        m_cluster_cube = cluster_cube;
        update_probability();
    }

    void set_number_pp(double number_pp){
        m_number_pp=number_pp;
    }

    std::complex<double> eval(const std::complex<double>& s) const override {
        return std::pow( expected_single_point(s),  m_number_pp );
    }

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    // L(s) = exp(N * eval_exponent(s)): the BPP counterpart of the PPP exponent, i.e. d log L(s) / dN.
    std::complex<double> eval_exponent(const std::complex<double>& s) const override {
        return std::log( expected_single_point(s) );
    }

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval_exponent(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};

} // namespace risgeo
