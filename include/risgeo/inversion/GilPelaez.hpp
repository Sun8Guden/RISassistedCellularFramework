#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "../integration/cuhre/Integration.hpp"
#include "../integration/cuhre/Region.hpp"

namespace risgeo {

// Inversion formulas for a real random variable X known only through its two-sided Laplace transform
// L(s) = E[exp(-s X)], evaluated on vertical lines Re(s) = a of the complex plane:
//
//   Gil-Pelaez               P(X > 0)               = 1/2      + 1/pi \int_0^\infty Im L(-iu)    / u du
//   exponentially tilted     E[exp(-a X) 1{X > 0}]  = L(a) / 2 + 1/pi \int_0^\infty Im L(a - iu) / u du
//
// The second line is Gil-Pelaez under the measure tilted by exp(-a X). X may take both signs, which is
// what happens to the "effective interference" theta (I + N) - R once reflected signals R are present.
//
// For a Rayleigh-faded direct link the coverage probability is P(H > X) with H ~ Exp(1), independent of X:
//   P(H > X) = E[exp(-X) 1{X > 0}] + P(X <= 0)
//            = 1/2 (1 + L(1)) + 1/pi \int_0^\infty Im[ L(1 - iu) - L(-iu) ] / u du,
// evaluated as a single integral so that the two oscillating terms cancel inside the integrand.
//
// The same formulas price options from a characteristic function (see examples/OptionPricingExample.cpp).
namespace GilPelaez {

    // (1/pi) \int_0^\infty Im[g(u)] / u du with double-exponential (exp-sinh) quadrature.
    template <typename F>
    double imaginary_part_integral(F g, double tolerance = 1e-6){
        boost::math::quadrature::exp_sinh<double> singular;
        auto integrand = [&](const double& u){
            std::complex<double> value = g(u);
            return value.imag() / (u * M_PI);
        };
        double error;
        double L1;
        return singular.integrate(integrand, tolerance, &error, &L1);
    }

    // P(X > 0)
    template <typename LT>
    double probability_positive(LT laplace, double tolerance = 1e-6){
        double cpv_part_value = imaginary_part_integral([&](double u){
            return laplace(std::complex<double>{0.0, - u});
        }, tolerance);
        return 0.5 + cpv_part_value;
    }

    // E[exp(-a X) 1{X > 0}] for a real a inside the strip where L converges.
    template <typename LT>
    double tilted_expectation_positive(LT laplace, double a, double tolerance = 1e-6){
        double cpv_part_value = imaginary_part_integral([&](double u){
            return laplace(std::complex<double>{a, - u});
        }, tolerance);
        return 0.5 * laplace(std::complex<double>{a, 0.0}).real() + cpv_part_value;
    }

    // P(H > X) with H ~ Exp(1) independent of X.
    template <typename LT>
    double exponential_exceeds(LT laplace, double tolerance = 1e-6){
        double cpv_part_value = imaginary_part_integral([&](double u){
            std::complex<double> cpv_part_one = laplace(std::complex<double>{1.0, - u});
            std::complex<double> cpv_part_two = laplace(std::complex<double>{0.0, - u});
            return cpv_part_one - cpv_part_two;
        }, tolerance);
        std::complex<double> real_part_compute = laplace(std::complex<double>{1.0, 0.0});
        double real_part_value = 0.5 * (1.0 + real_part_compute.real());
        return cpv_part_value + real_part_value;
    }

    // \int_0^{theta_max} P(H > X_theta) / (1 + theta) d theta for a family X_theta, e.g. the ergodic rate
    // E[ln(1 + SINR)] when X_theta is the effective interference at SINR threshold theta.
    // weighted_laplace(s, theta) must return L_theta(s) / (1 + theta). The Gil-Pelaez integral is taken in
    // log z (z = exp(u), u in [-log_z_width/2, log_z_width/2]), which turns dz / z into du and removes the
    // singularity at z = 0, so (theta, u) is a smooth 2D integral for the adaptive cubature.
    template <typename F>
    double integrated_exponential_exceeds(F weighted_laplace, double theta_max, double log_z_width,
            IntegrationStrategy strategy = IntegrationStrategy::parallel){
        std::array<double, 2> center2d{theta_max/2.0,  0.0};
        std::array<double, 2> width2d{theta_max,  log_z_width};
        Region<2> cur_region2d = Region<2>(center2d, width2d);

        double cpv_part_value{ 0.0 }, real_part_value{ 0.0 };
        double estimated_error_part_1{ 0.0 }, estimated_error_part_2{ 0.0 };
        int number_f_eval_part1{ 0 };

        auto f_cpv = [&](const std::array<double, 2>& input2d) {
            double z = - std::exp(input2d.at(1));
            std::complex<double> cpv_part_one = weighted_laplace( std::complex<double>{1.0, z}, input2d.at(0) );
            std::complex<double> cpv_part_two = weighted_laplace( std::complex<double>{0.0, z}, input2d.at(0) );
            std::complex<double> cpv =  cpv_part_one - cpv_part_two;
            double cpv_val = cpv.imag();
            double ret_val = cpv_val / M_PI;
            return ret_val;
        };
        cpv_part_value = Integration<2>::integrate(f_cpv, cur_region2d, estimated_error_part_1, \
                number_f_eval_part1, strategy, 1e-6, 500, 50000);

        auto f_real = [&](const double threshold ){
            std::complex<double> real_part_compute =  weighted_laplace(std::complex<double>{1.0, 0.0}, threshold);
            double real_part = 0.5 * (1.0/(1.0+threshold) + real_part_compute.real());
            return real_part;
        };
        real_part_value =  boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f_real, 0.0, theta_max, 10, 1e-6, &estimated_error_part_2);

        return real_part_value + cpv_part_value;
    }

} // namespace GilPelaez

} // namespace risgeo
