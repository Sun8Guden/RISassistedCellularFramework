#pragma once
#include <complex>
#include <cmath>
#include "Parameter.hpp"
#include "CoverageProbability.hpp"
#include "../laplace/LaplaceTransform.hpp"
#include "../inversion/GilPelaez.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace risgeo {

// Ergodic rate E[ln(1 + SINR)] = \int_0^\infty P(SINR > theta) / (1 + theta) d theta, in nat/s/Hz.
// With reflected signals the Gil-Pelaez term turns this into a 2D integral over (theta, log z),
// see GilPelaez::integrated_exponential_exceeds.

class ErgodicRate {

    public:

    static double calculate_rayleigh(const double distance, const LaplaceTransform& interference, const double& noise){
        auto coverage_probability = [&](double threshold){
            double coverage_val =  CoverageProbability::calculate_rayleigh(distance, threshold, interference, noise);
            return coverage_val /(1.0 + threshold);
        };
        double estimated_error;
        double estimated_ergodic_rate = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(coverage_probability, 0.0, \
                Parameter::infinity_large_threshold, 10, 1e-6, &estimated_error);
        return estimated_ergodic_rate;
    }

    static double calculate_rayleigh_positive(const double distance, const LaplaceTransform& interference, \
            const double noise, const LaplaceTransform& reflected_signals){

        auto f_laplace_eval = [&](const std::complex<double>& s, const double threshold){
            std::complex<double> coverage_probability = interference.eval( s * threshold/PathLoss::direct_nlos(distance) ) *\
            reflected_signals.eval( s / PathLoss::direct_nlos(distance)) * \
            GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
            return coverage_probability / (1.0 + threshold);
        };

        return GilPelaez::integrated_exponential_exceeds(f_laplace_eval, Parameter::infinity_large_threshold,
                Parameter::infinity_cauchy_principal_value);
    }
};

} // namespace risgeo
