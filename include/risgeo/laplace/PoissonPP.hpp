#pragma once
#include <complex>
#include <limits>
#include <functional>
#include <cmath>
#include "LaplaceTransform.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "../integration/genz_malik/Cube.hpp"
#include "../integration/genz_malik/GM2D.hpp"

namespace risgeo {

// We introduce two classes for 2D Poisson point processes: one for centered scenario
// where we can reduce one dimension whereas the other for non-centered scenario.
// They are mainly of interest for stochastic geometry research.

// For a PPP of density lambda, L(s) = exp(lambda * eval_exponent(s)), so eval_exponent(s)
// is also d log L(s) / d lambda, which is what the density derivatives in the research code use.


class PoissonPP2D_Reduced: public LaplaceTransform {
    double m_density_pp;
    std::function<std::complex<double>(std::complex<double>, double)> mf_laplace_function_1D;
    double m_region_left;
    double m_region_right;

    std::complex<double> radial_integral(const std::complex<double>& s) const {
        double estimated_error;
        auto f_integrand = [&](const double distance){
            return distance * (1.0 - mf_laplace_function_1D(s, distance));
        };
        return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f_integrand, \
        m_region_left,  m_region_right, 10, 1e-6, &estimated_error);
    }

    public:
    PoissonPP2D_Reduced(const std::function<std::complex<double>(std::complex<double>, double)>& \
        laplace_function_1D, double density_pp, double x_left, double x_right = std::numeric_limits<double>::infinity()\
        ):m_density_pp(density_pp), mf_laplace_function_1D(laplace_function_1D), m_region_left(x_left), m_region_right(x_right){};

    void set_functor(const std::function<std::complex<double>(std::complex<double>, double)>& laplace_function_1D){
        mf_laplace_function_1D = laplace_function_1D;
    };

    void set_region(double x_left, double x_right = std::numeric_limits<double>::infinity() ){
        m_region_left = x_left;
        m_region_right = x_right;
    };

    std::complex<double> eval(const std::complex<double>& s) const override {
        std::complex<double> integral_value = radial_integral(s);
        return std::exp( - 2.0 * M_PI * integral_value * m_density_pp);
    };

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    std::complex<double> eval_exponent(const std::complex<double>& s) const override {
        std::complex<double> integral_value = radial_integral(s);
        return  - 2.0 * M_PI * integral_value;
    };

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval_exponent(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};


class PoissonPP2D : public LaplaceTransform {
        std::function<std::complex<double>(std::complex<double>, double, double)> mf_laplace_transform_2D;
        CUBE::Cube<double, 2> m_cluster_cube;
        double m_density_pp;

        std::complex<double> cluster_integral(const std::complex<double>& s) const {
            double estimated_error{0.0};
            unsigned int number_evaluation{0};
            auto f_integrand = [&](const double y, const double theta){
                return y * (1.0 - mf_laplace_transform_2D(s, y, theta));
            };
            return GM::GM2D<double>::integrate(f_integrand, \
                m_cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        }

    public:
        PoissonPP2D(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D, \
            double density_pp, const CUBE::Cube<double, 2>& cluster_cube): mf_laplace_transform_2D(laplace_function_2D),\
            m_cluster_cube(cluster_cube), m_density_pp(density_pp){}

        void set_functor(const std::function<std::complex<double>(std::complex<double>, double, double)>& laplace_function_2D){
            mf_laplace_transform_2D = laplace_function_2D;
        }

        void set_region(const CUBE::Cube<double, 2>& cluster_cube){
            m_cluster_cube = cluster_cube;
        }

    std::complex<double> eval(const std::complex<double>& s) const override {
        std::complex<double> integral_value = cluster_integral(s);
        return std::exp( - integral_value * m_density_pp );
    }

    double eval(const double& s) const override {
        std::complex<double> complex_wrapper = eval(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

    std::complex<double> eval_exponent(const std::complex<double>& s) const override {
        std::complex<double> integral_value = cluster_integral(s);
        return  - integral_value;
    }

    double eval_exponent(const double& s) const override {
        std::complex<double> complex_wrapper = eval_exponent(std::complex<double>{s, 0.0});
        return complex_wrapper.real();
    };

};

} // namespace risgeo
