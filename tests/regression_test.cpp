// Regression and consistency checks for the framework (run by `make test` and `ctest`).
// Reference values are the results documented in examples/Readme.md and in the literature cited below.
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include "../include/risgeo/laplace/PoissonPP.hpp"
#include "../include/risgeo/laplace/BinomialPP.hpp"
#include "../include/risgeo/laplace/Exponential.hpp"
#include "../include/risgeo/laplace/ChiSquare.hpp"
#include "../include/risgeo/geometry/CoverageProbability.hpp"
#include "../include/risgeo/geometry/ErgodicRate.hpp"
#include "../include/risgeo/integration/cuhre/Integration.hpp"
#include "../include/risgeo/integration/genz_malik/GM2D.hpp"
#include "../include/risgeo/utils/pathloss.hpp"
#include "../examples/finance/pricing.hpp"

using namespace risgeo;

namespace {

int g_failures = 0;

void check(const std::string& name, double value, double expected, double rel_tol){
    double rel_error = std::abs(value - expected) / std::abs(expected);
    bool ok = rel_error <= rel_tol;
    if (!ok) g_failures++;
    std::printf("[%s] %-52s %.10g (expected %.10g, rel. error %.1e)\n", ok ? "PASS" : "FAIL", name.c_str(), value, expected, rel_error);
}

// Scenario of the examples: BS density 10/km^2, typical link of 100 m, RIS ring of [10, 25] m with 5 RISs.
constexpr double noise_power = 1e-13;
constexpr double BS_density = 1e-5;
constexpr double distance_typical_BS = 100.0;
constexpr double RIS_elements = 200.0;
constexpr double radius_min = 10.0;
constexpr double radius_max = 25.0;

std::complex<double> lt_interference_at_r(std::complex<double> s, double r){
    return Exponential::eval(s, PathLoss::direct_nlos(r));
}

std::complex<double> lt_reflected_signal_at_y_theta(std::complex<double> s, double y, double theta){
    double chi_variance_parameter = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance;
    double non_center_parameter = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean / Parameter::zeta_variance;
    chi_variance_parameter *= PathLoss::direct_los(y) * PathLoss::reflected_los(distance_typical_BS, y, theta);
    return ChiSquare::eval_non_central(-s, non_center_parameter, chi_variance_parameter);
}

void test_examples(){
    auto cluster_cube = CUBE::make_cube_2D<double>(radius_min, radius_max, 0.0, 2.0 * M_PI);
    double RIS_density = 5.0 / (M_PI * (radius_max * radius_max - radius_min * radius_min));
    PoissonPP2D_Reduced interference(lt_interference_at_r, BS_density, distance_typical_BS);
    PoissonPP2D reflected_signals(lt_reflected_signal_at_y_theta, RIS_density, cluster_cube);

    check("coverage probability", CoverageProbability::calculate_rayleigh(distance_typical_BS, 1.0, interference, noise_power),
          0.660685, 1e-5);
    check("coverage probability with RISs", CoverageProbability::calculate_rayleigh_positive(distance_typical_BS, 1.0, noise_power,
          interference, reflected_signals), 0.777182, 1e-5);
    check("ergodic rate", ErgodicRate::calculate_rayleigh(distance_typical_BS, interference, noise_power), 1.11374, 1e-5);
    // The adaptive 2D integration subdivides differently with and without OpenMP, hence the looser tolerance.
    check("ergodic rate with RISs", ErgodicRate::calculate_rayleigh_positive(distance_typical_BS, interference, noise_power,
          reflected_signals), 1.27658, 1e-3);
}

void test_integrators(){
    // \int_0^1 \int_0^1 exp(-x - y) dx dy = (1 - 1/e)^2
    const double exact = (1.0 - std::exp(-1.0)) * (1.0 - std::exp(-1.0));

    Region<2> region({0.5, 0.5}, {1.0, 1.0});
    double estimated_error{0.0};
    int num_func_eval{0};
    double cuhre = Integration<2>::integrate([](const std::array<double, 2>& x){ return std::exp(-x[0] - x[1]); },
                                             region, estimated_error, num_func_eval);
    check("Cuhre 2D integral of exp(-x-y)", cuhre, exact, 1e-9);

    unsigned int num_f_eval{0};
    double gm_error{0.0};
    double genz_malik = GM::GM2D<double>::integrate([](double x, double y){ return std::exp(-x - y); },
                                                    CUBE::make_cube_2D<double>(0.0, 1.0, 0.0, 1.0), 1e-9, gm_error, 5000, num_f_eval);
    check("Genz-Malik 2D integral of exp(-x-y)", genz_malik, exact, 1e-9);
}

void test_laplace_derivatives(){
    // The *_D1 functions must match a central finite difference of the transform.
    const std::complex<double> s{0.3, -0.7};
    const double h = 1e-6, variance = 0.8, non_central = 2.5;
    auto finite_difference = [&](auto f){ return std::abs((f(s + h) - f(s - h)) / (2.0 * h)); };

    auto exp_fd = finite_difference([&](std::complex<double> x){ return Exponential::eval(x, variance); });
    check("Exponential::eval_D1 vs finite difference", std::abs(Exponential::eval_D1(s, variance)), exp_fd, 1e-6);
    auto chi_fd = finite_difference([&](std::complex<double> x){ return ChiSquare::eval_central(x, variance); });
    check("ChiSquare::eval_central_D1 vs finite difference", std::abs(ChiSquare::eval_central_D1(s, variance)), chi_fd, 1e-6);
    auto nc_fd = finite_difference([&](std::complex<double> x){ return ChiSquare::eval_non_central(x, non_central, variance); });
    check("ChiSquare::eval_non_central_D1 vs finite difference", std::abs(ChiSquare::eval_non_central_D1(s, non_central, variance)), nc_fd, 1e-6);

    // Sign convention: dL/ds is negative for a positive random variable at real s > 0.
    check("sign of ChiSquare::eval_central_D1", ChiSquare::eval_central_D1(0.5, variance) < 0.0 ? 1.0 : 0.0, 1.0, 0.0);
}

void test_point_processes(){
    // L(s) = exp(intensity * eval_exponent(s)) for both overloads.
    PoissonPP2D_Reduced interference(lt_interference_at_r, BS_density, distance_typical_BS);
    const double s_real = 1.0 / PathLoss::direct_nlos(distance_typical_BS);
    check("PoissonPP2D_Reduced exponent (real s)", std::exp(BS_density * interference.eval_exponent(s_real)),
          interference.eval(s_real), 1e-12);
    const std::complex<double> s_complex{s_real, -0.5 * s_real};
    check("PoissonPP2D_Reduced exponent (complex s)", std::abs(std::exp(BS_density * interference.eval_exponent(s_complex))),
          std::abs(interference.eval(s_complex)), 1e-12);

    auto cluster_cube = CUBE::make_cube_2D<double>(radius_min, radius_max, 0.0, 2.0 * M_PI);
    const double number_RIS = 5.0;
    BinomialPP reflected_BPP(lt_reflected_signal_at_y_theta, number_RIS, cluster_cube);
    const std::complex<double> s_signal{0.0, -0.2 * s_real};
    check("BinomialPP exponent", std::abs(std::exp(number_RIS * reflected_BPP.eval_exponent(s_signal))),
          std::abs(reflected_BPP.eval(s_signal)), 1e-12);

    // A point with a trivial transform (L = 1) must give L = 1 for the whole process: the ring density is normalised.
    BinomialPP trivial_BPP([](std::complex<double>, double, double){ return std::complex<double>{1.0, 0.0}; }, number_RIS, cluster_cube);
    check("BinomialPP normalisation", trivial_BPP.eval(0.1), 1.0, 1e-9);
}

void test_option_pricing(){
    // The inversion engine on problems with closed-form or published answers.
    const finance::Market market{100.0, 0.03, 1.0};
    const finance::BlackScholes bs{0.2};
    for (double strike : {80.0, 100.0, 120.0}) {
        check("Black-Scholes call K=" + std::to_string(int(strike)) + " vs closed form", finance::call(bs, market, strike),
              finance::black_scholes_call(market, strike, 0.2), 1e-12);
        check("Black-Scholes digital K=" + std::to_string(int(strike)) + " vs closed form", finance::digital_call(bs, market, strike),
              finance::black_scholes_digital_call(market, strike, 0.2), 1e-12);
    }

    // Fang and Oosterlee, SIAM J. Sci. Comput. 31 (2008): Heston reference price 5.785155434.
    const finance::Heston benchmark{1.5768, 0.0398, 0.5751, -0.5711, 0.0175};
    check("Heston call vs Fang & Oosterlee (2008)", finance::call(benchmark, finance::Market{100.0, 0.0, 1.0}, 100.0),
          5.785155434, 1e-9);

    // Martingale property of the characteristic function: E[S_T] = S_0 e^{rT}.
    const finance::Heston heston{1.5, 0.04, 0.3, -0.7, 0.04};
    check("Heston forward phi(-i) = S0 exp(rT)", heston.log_price_cf({0.0, -1.0}, market).real(),
          market.spot * std::exp(market.rate * market.maturity), 1e-12);

    // Vanishing vol-of-vol with v0 = theta reduces Heston to Black-Scholes with sigma = sqrt(theta).
    check("Heston with xi -> 0 vs Black-Scholes", finance::call(finance::Heston{1.5, 0.04, 1e-5, -0.7, 0.04}, market, 110.0),
          finance::black_scholes_call(market, 110.0, 0.2), 1e-4);
}

} // namespace

int main(){
    test_integrators();
    test_laplace_derivatives();
    test_point_processes();
    test_option_pricing();
    test_examples();
    std::printf("\n%s: %d failure(s)\n", g_failures == 0 ? "OK" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
