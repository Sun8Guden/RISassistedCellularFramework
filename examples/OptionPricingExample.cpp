// The inversion engine behind the network metrics, applied to option pricing.
//
// 1. Black-Scholes: Gil-Pelaez prices against the closed form.
// 2. Heston (stochastic volatility, no closed form): the published reference price of Fang & Oosterlee (2008),
//    then an implied-volatility smile cross-checked against an independent Monte Carlo simulation.
//
// Usage: OptionPricingExample [output_dir]   (writes heston_smile.csv when a directory is given)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "finance/pricing.hpp"

namespace {

double seconds_since(std::chrono::steady_clock::time_point start){
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// Terminal prices of a full-truncation Euler scheme for Heston, with antithetic variates.
// Paths are simulated in chunks with their own seeds, so the result does not depend on the number of threads.
std::vector<double> heston_terminal_prices(const finance::Heston& h, const finance::Market& m, long pairs, int steps){
    const long chunk = 5000;
    const long chunks = (pairs + chunk - 1) / chunk;
    std::vector<double> terminal(2 * pairs);
    const double dt = m.maturity / steps;
    #pragma omp parallel for schedule(dynamic)
    for (long c = 0; c < chunks; c++) {
        std::mt19937_64 rng(7 + c);
        std::normal_distribution<double> gaussian(0.0, 1.0);
        for (long p = c * chunk; p < std::min(pairs, (c + 1) * chunk); p++) {
            double log_s[2] = {std::log(m.spot), std::log(m.spot)};
            double v[2] = {h.v0, h.v0};
            for (int k = 0; k < steps; k++) {
                const double z_v = gaussian(rng);
                const double z_s = h.rho * z_v + std::sqrt(1.0 - h.rho * h.rho) * gaussian(rng);
                for (int a = 0; a < 2; a++) {
                    const double sign = a == 0 ? 1.0 : -1.0;
                    const double v_plus = std::max(v[a], 0.0);
                    log_s[a] += (m.rate - 0.5 * v_plus) * dt + std::sqrt(v_plus * dt) * sign * z_s;
                    v[a] += h.kappa * (h.theta - v_plus) * dt + h.xi * std::sqrt(v_plus * dt) * sign * z_v;
                }
            }
            terminal[2 * p] = std::exp(log_s[0]);
            terminal[2 * p + 1] = std::exp(log_s[1]);
        }
    }
    return terminal;
}

} // namespace

int main(int argc, char** argv){
    const std::string output_dir = argc > 1 ? argv[1] : "";
    const finance::Market market{100.0, 0.03, 1.0};
    const std::vector<double> strikes{70, 80, 90, 95, 100, 105, 110, 120, 130};

    // ---------------- Black-Scholes ----------------
    const finance::BlackScholes bs{0.2};
    std::printf("Black-Scholes (S0 = 100, r = 3%%, sigma = 20%%, T = 1y): Gil-Pelaez vs. closed form\n");
    std::printf("%8s %14s %10s %14s %10s\n", "strike", "call", "abs. error", "digital", "abs. error");
    double max_error = 0.0;
    auto start = std::chrono::steady_clock::now();
    for (double k : strikes) {
        const double call = finance::call(bs, market, k), digital = finance::digital_call(bs, market, k);
        const double e_call = std::abs(call - finance::black_scholes_call(market, k, bs.volatility));
        const double e_digital = std::abs(digital - finance::black_scholes_digital_call(market, k, bs.volatility));
        max_error = std::max({max_error, e_call, e_digital});
        std::printf("%8.1f %14.8f %10.1e %14.8f %10.1e\n", k, call, e_call, digital, e_digital);
    }
    std::printf("Largest error %.1e, %.0f microseconds per price (call + digital)\n\n",
                max_error, 1e6 * seconds_since(start) / strikes.size());

    // ---------------- Heston: published reference ----------------
    // F. Fang and C. W. Oosterlee, A novel pricing method for European options based on Fourier-cosine
    // series expansions, SIAM J. Sci. Comput. 31 (2008): reference call price 5.785155434.
    const double reference = 5.785155434;
    const double benchmark = finance::call(finance::Heston{1.5768, 0.0398, 0.5751, -0.5711, 0.0175}, finance::Market{100.0, 0.0, 1.0}, 100.0);
    std::printf("Heston benchmark of Fang & Oosterlee (2008): %.10f (published %.9f, difference %.1e)\n\n",
                benchmark, reference, std::abs(benchmark - reference));

    // ---------------- Heston: smile and Monte Carlo cross-check ----------------
    const finance::Heston heston{1.5, 0.04, 0.3, -0.7, 0.04};
    std::printf("Heston (kappa = 1.5, theta = 0.04, xi = 0.3, rho = -0.7, v0 = 0.04): Gil-Pelaez vs. Monte Carlo\n");
    start = std::chrono::steady_clock::now();
    std::vector<double> prices;
    for (double k : strikes) prices.push_back(finance::call(heston, market, k));
    const double inversion_seconds = seconds_since(start);

    const long pairs = 400000;
    const int steps = 1000;
    start = std::chrono::steady_clock::now();
    const std::vector<double> terminal = heston_terminal_prices(heston, market, pairs, steps);
    const double mc_seconds = seconds_since(start);

    std::ofstream smile;
    if (!output_dir.empty()) {
        smile.open(output_dir + "/heston_smile.csv");
        smile << "strike,inversion_price,inversion_implied_vol,mc_price,mc_ci_low,mc_ci_high,mc_implied_vol_low,mc_implied_vol_high\n";
    }
    std::printf("%8s %12s %12s %20s %8s %10s\n", "strike", "inversion", "Monte Carlo", "95% CI", "z", "impl. vol");
    const double discount = std::exp(-market.rate * market.maturity);
    double max_z = 0.0;
    for (std::size_t j = 0; j < strikes.size(); j++) {
        const double k = strikes[j];
        // Average each antithetic pair first so that the standard error accounts for their correlation.
        double sum = 0.0, sum_sq = 0.0;
        for (long p = 0; p < pairs; p++) {
            const double payoff = 0.5 * discount * (std::max(terminal[2 * p] - k, 0.0) + std::max(terminal[2 * p + 1] - k, 0.0));
            sum += payoff;
            sum_sq += payoff * payoff;
        }
        const double mean = sum / pairs;
        const double standard_error = std::sqrt((sum_sq / pairs - mean * mean) / pairs);
        const double z = (mean - prices[j]) / standard_error;
        max_z = std::max(max_z, std::abs(z));
        const double implied = finance::implied_volatility(market, k, prices[j]);
        std::printf("%8.1f %12.6f %12.6f [%8.4f, %8.4f] %8.2f %9.2f%%\n", k, prices[j], mean,
                    mean - 1.96 * standard_error, mean + 1.96 * standard_error, z, 100.0 * implied);
        if (smile.is_open()) {
            smile << k << "," << prices[j] << "," << implied << "," << mean << "," << mean - 1.96 * standard_error << ","
                  << mean + 1.96 * standard_error << "," << finance::implied_volatility(market, k, mean - 1.96 * standard_error) << ","
                  << finance::implied_volatility(market, k, mean + 1.96 * standard_error) << "\n";
        }
    }
    std::printf("Inversion: %.1f ms for %zu strikes. Monte Carlo: %.1f s for %ld paths x %d steps (%ld antithetic pairs).\n",
                1e3 * inversion_seconds, strikes.size(), mc_seconds, 2 * pairs, steps, pairs);
    std::printf("Largest |MC - inversion| / standard error: %.2f (the Euler scheme also has a small discretisation bias)\n", max_z);
    return max_error < 1e-7 && std::abs(benchmark - reference) < 1e-8 ? 0 : 1;
}
