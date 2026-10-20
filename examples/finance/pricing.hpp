#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include "../../include/risgeo/inversion/GilPelaez.hpp"

// European option pricing with the same Gil-Pelaez inversion that computes the network metrics.
//
// With X = ln(S_T / K) under the risk-neutral measure and L(s) = E[exp(-s X)] = K^s phi(i s), where
// phi(u) = E[exp(i u ln S_T)] is the characteristic function of the log-price:
//   digital call = e^{-rT} P(X > 0)                                   GilPelaez::probability_positive
//   call         = e^{-rT} K ( E[e^X 1{X > 0}] - P(X > 0) )           GilPelaez::tilted_expectation_positive, a = -1
// The tilt a = -1 is the change to the stock-price numeraire, exactly as a = 1 is the Rayleigh-fading tilt
// of the coverage probability.
namespace finance {

struct Market {
    double spot;
    double rate;
    double maturity;
};

// ln S_T ~ N(ln S_0 + (r - sigma^2 / 2) T, sigma^2 T)
struct BlackScholes {
    double volatility;

    std::complex<double> log_price_cf(std::complex<double> u, const Market& m) const {
        const std::complex<double> i{0.0, 1.0};
        const double variance = volatility * volatility * m.maturity;
        const double mean = std::log(m.spot) + m.rate * m.maturity - 0.5 * variance;
        return std::exp(i * u * mean - 0.5 * variance * u * u);
    }
};

// dS = r S dt + sqrt(v) S dW_1,  dv = kappa (theta - v) dt + xi sqrt(v) dW_2,  d<W_1, W_2> = rho dt.
// Characteristic function in the "little Heston trap" form of Albrecher et al. (2007), which keeps the
// complex logarithm on its principal branch.
struct Heston {
    double kappa, theta, xi, rho, v0;

    std::complex<double> log_price_cf(std::complex<double> u, const Market& m) const {
        const std::complex<double> i{0.0, 1.0};
        const double T = m.maturity;
        const std::complex<double> beta = kappa - rho * xi * i * u;
        const std::complex<double> d = std::sqrt(beta * beta + xi * xi * (i * u + u * u));
        const std::complex<double> g = (beta - d) / (beta + d);
        const std::complex<double> exp_dT = std::exp(-d * T);
        const std::complex<double> C = kappa * theta / (xi * xi) * ((beta - d) * T - 2.0 * std::log((1.0 - g * exp_dT) / (1.0 - g)));
        const std::complex<double> D = (beta - d) / (xi * xi) * (1.0 - exp_dT) / (1.0 - g * exp_dT);
        return std::exp(i * u * (std::log(m.spot) + m.rate * T) + C + D * v0);
    }
};

// L(s) = E[exp(-s ln(S_T / K))]
template <typename Model>
auto log_moneyness_laplace(const Model& model, const Market& market, double strike){
    return [=](std::complex<double> s){
        const std::complex<double> i{0.0, 1.0};
        return std::exp(s * std::log(strike)) * model.log_price_cf(i * s, market);
    };
}

template <typename Model>
double digital_call(const Model& model, const Market& market, double strike, double tolerance = 1e-10){
    auto laplace = log_moneyness_laplace(model, market, strike);
    return std::exp(-market.rate * market.maturity) * risgeo::GilPelaez::probability_positive(laplace, tolerance);
}

template <typename Model>
double call(const Model& model, const Market& market, double strike, double tolerance = 1e-10){
    auto laplace = log_moneyness_laplace(model, market, strike);
    const double in_the_money = risgeo::GilPelaez::probability_positive(laplace, tolerance);
    const double share_weighted = risgeo::GilPelaez::tilted_expectation_positive(laplace, -1.0, tolerance);
    return std::exp(-market.rate * market.maturity) * strike * (share_weighted - in_the_money);
}

template <typename Model>
double put(const Model& model, const Market& market, double strike, double tolerance = 1e-10){
    return call(model, market, strike, tolerance) - market.spot + strike * std::exp(-market.rate * market.maturity);
}

// ---------------- Closed forms and implied volatility ----------------

inline double normal_cdf(double x){ return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

inline double black_scholes_call(const Market& m, double strike, double volatility){
    const double sqrt_T = std::sqrt(m.maturity);
    const double d1 = (std::log(m.spot / strike) + (m.rate + 0.5 * volatility * volatility) * m.maturity) / (volatility * sqrt_T);
    const double d2 = d1 - volatility * sqrt_T;
    return m.spot * normal_cdf(d1) - strike * std::exp(-m.rate * m.maturity) * normal_cdf(d2);
}

inline double black_scholes_digital_call(const Market& m, double strike, double volatility){
    const double sqrt_T = std::sqrt(m.maturity);
    const double d2 = (std::log(m.spot / strike) + (m.rate - 0.5 * volatility * volatility) * m.maturity) / (volatility * sqrt_T);
    return std::exp(-m.rate * m.maturity) * normal_cdf(d2);
}

// Black-Scholes implied volatility by bisection (the call price is increasing in the volatility).
inline double implied_volatility(const Market& m, double strike, double price){
    double low = 1e-4, high = 3.0;
    for (int iteration = 0; iteration < 200 && high - low > 1e-12; iteration++) {
        const double mid = 0.5 * (low + high);
        (black_scholes_call(m, strike, mid) < price ? low : high) = mid;
    }
    return 0.5 * (low + high);
}

} // namespace finance
