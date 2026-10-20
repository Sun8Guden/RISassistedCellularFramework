#pragma once
#define _USE_MATH_DEFINES // for using the constant PI in C++
#include <cmath>
#include <complex>
#include <functional>
#include "Parameter.hpp"
#include "../laplace/LaplaceTransform.hpp"
#include "../laplace/PoissonPP.hpp"
#include "../laplace/GaussianNoise.hpp"
#include "../inversion/GilPelaez.hpp"
#include "../utils/pathloss.hpp"

namespace risgeo {

// We want to provide a general methodology for computing coverage probability for multiple common direct links,
// with or without the reflected signals.
//
// The typical user is served by a BS at distance r over a Rayleigh-faded link, so it is covered when
// H > X with H ~ Exp(1) and X = theta (I + N) / l(r) - R / l(r), where I is the interference, N the noise and
// R the power reflected by RISs. Without reflected signals X >= 0 and P(cov) = L_X(1) in closed form.
// With reflected signals X can be negative and the coverage probability needs the tilted Gil-Pelaez
// inversion of GilPelaez::exponential_exceeds:
//   P(cov) = 1/2 (1 + Re L_X(1)) + 1/pi \int_0^\infty Im[ L_X(1 - iz) - L_X(-iz) ] / z dz.
namespace CoverageProbability {

    inline double calculate_rayleigh(const double distance, const double threshold, \
            const LaplaceTransform& interference, const double noise){

        double coverage_probability = interference.eval( threshold  / PathLoss::direct_nlos(distance) )*\
            GaussianNoise::eval(threshold/(PathLoss::direct_nlos(distance) * Parameter::antenna_gain*Parameter::transmission_power), noise);

        return coverage_probability;

    }

    inline double calculate_rayleigh_positive(const double distance, const double threshold, const double noise, \
            const LaplaceTransform& interference, const LaplaceTransform& reflected_signals)
    {
        auto f_laplace_eval = [&](const std::complex<double>& s){
            std::complex<double> coverage_probability = interference.eval( s * threshold/PathLoss::direct_nlos(distance) )*\
            reflected_signals.eval( s / PathLoss::direct_nlos(distance) ) *\
            GaussianNoise::eval(s * threshold/(PathLoss::direct_nlos(distance)*Parameter::antenna_gain*Parameter::transmission_power), noise);
            return coverage_probability;
        };

        return GilPelaez::exponential_exceeds(f_laplace_eval, Parameter::expected_error);
    }

    // Same as calculate_rayleigh_positive for a user-supplied L_X(s; distance, threshold).
    inline double calculate_laplace_positive(std::function<std::complex<double>(std::complex<double>, double, double)> f_laplace_eval, const double distance, const double threshold){
        return GilPelaez::exponential_exceeds([&](const std::complex<double>& s){
            return f_laplace_eval(s, distance, threshold);
        }, Parameter::expected_error);
    }
}

} // namespace risgeo
