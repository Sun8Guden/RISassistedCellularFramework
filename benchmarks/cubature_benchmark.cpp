// Benchmarks the adaptive cubature (risgeo::Integration<N>).
//
// 1. Overhead-dominated: a cheap 3D Genz "Gaussian" integrand with a known integral, run to a fixed
//    budget of function evaluations, so the cost of managing the subregions is what is measured.
// 2. Integrand-dominated: the ergodic rate with RIS-reflected signals (examples/ErgodicRateRISExample.cpp).
// 3. Convergence: relative error vs. number of evaluations for the cubature and for plain Monte Carlo,
//    written to a CSV file when an output path is given (used by scripts/make_figures.py).
//
// Usage: cubature_benchmark [convergence.csv]
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <random>
#include "../include/risgeo/risgeo.hpp"

using namespace risgeo;

namespace {

double seconds_since(std::chrono::steady_clock::time_point start){
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// Genz Gaussian test function on [0, 1]^N: exp(-sum_i a_i^2 (x_i - w_i)^2).
template <int N>
struct GenzGaussian {
    std::array<double, N> a;
    std::array<double, N> w;

    double operator()(const std::array<double, N>& x) const {
        double exponent = 0.0;
        for (int i = 0; i < N; i++) exponent += a[i] * a[i] * (x[i] - w[i]) * (x[i] - w[i]);
        return std::exp(-exponent);
    }

    double exact() const {
        double value = 1.0;
        for (int i = 0; i < N; i++) value *= std::sqrt(M_PI) / (2.0 * a[i]) * (std::erf(a[i] * (1.0 - w[i])) + std::erf(a[i] * w[i]));
        return value;
    }
};

template <int N>
Region<N> unit_cube(){
    std::array<double, N> center, width;
    center.fill(0.5);
    width.fill(1.0);
    return Region<N>(center, width);
}

void overhead_benchmark(){
    GenzGaussian<3> f{{5.0, 10.0, 20.0}, {0.3, 0.5, 0.7}};
    std::printf("Overhead-dominated: 3D Genz Gaussian, fixed evaluation budget\n");
    std::printf("%12s %12s %14s %14s\n", "evaluations", "time [s]", "ns / eval", "rel. error");
    for (int budget : {10000, 100000, 1000000}) {
        Region<3> region = unit_cube<3>();
        double estimated_error{0.0};
        int num_func_eval{0};
        auto start = std::chrono::steady_clock::now();
        // A tolerance of 0 forces the full budget.
        double value = Integration<3>::integrate(f, region, estimated_error, num_func_eval, IntegrationStrategy::serial, 0.0, 0, budget);
        double elapsed = seconds_since(start);
        std::printf("%12d %12.4f %14.1f %14.2e\n", num_func_eval, elapsed, 1e9 * elapsed / num_func_eval,
                    std::abs(value - f.exact()) / f.exact());
    }
}

void workload_benchmark(){
    const double noise_power = 1e-13, BS_density = 1e-5, distance = 100.0;
    const double RIS_elements = 200.0, radius_min = 10.0, radius_max = 25.0;
    auto cluster_cube = CUBE::make_cube_2D<double>(radius_min, radius_max, 0.0, 2.0 * M_PI);
    const double RIS_density = 5.0 / (M_PI * (radius_max * radius_max - radius_min * radius_min));

    PoissonPP2D_Reduced interference([](std::complex<double> s, double r){
        return Exponential::eval(s, PathLoss::direct_nlos(r));
    }, BS_density, distance);
    PoissonPP2D reflected_signals([&](std::complex<double> s, double y, double theta){
        double variance = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance
                          * PathLoss::direct_los(y) * PathLoss::reflected_los(distance, y, theta);
        double non_centrality = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean / Parameter::zeta_variance;
        return ChiSquare::eval_non_central(-s, non_centrality, variance);
    }, RIS_density, cluster_cube);

    auto start = std::chrono::steady_clock::now();
    double rate = ErgodicRate::calculate_rayleigh_positive(distance, interference, noise_power, reflected_signals);
    std::printf("\nIntegrand-dominated: ergodic rate with RISs = %.6f nat/s/Hz in %.3f s\n", rate, seconds_since(start));
}

void convergence_study(const char* path){
    GenzGaussian<2> f{{8.0, 8.0}, {0.4, 0.6}};
    const double exact = f.exact();
    std::ofstream out(path);
    out << "method,evaluations,relative_error\n";
    for (int budget = 200; budget <= 200000; budget = int(budget * 1.5)) {
        Region<2> region = unit_cube<2>();
        double estimated_error{0.0};
        int num_func_eval{0};
        double value = Integration<2>::integrate(f, region, estimated_error, num_func_eval, IntegrationStrategy::serial, 0.0, 0, budget);
        out << "cubature," << num_func_eval << "," << std::abs(value - exact) / exact << "\n";

        // Plain Monte Carlo with the same budget, averaged over repetitions for a stable curve.
        std::mt19937_64 rng(budget);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double mean_abs_error = 0.0;
        const int repetitions = 20;
        for (int r = 0; r < repetitions; r++) {
            double sum = 0.0;
            for (int i = 0; i < budget; i++) sum += f({uniform(rng), uniform(rng)});
            mean_abs_error += std::abs(sum / budget - exact) / exact / repetitions;
        }
        out << "monte_carlo," << budget << "," << mean_abs_error << "\n";
    }
    std::printf("\nConvergence data written to %s\n", path);
}

} // namespace

int main(int argc, char** argv){
    overhead_benchmark();
    workload_benchmark();
    if (argc > 1) convergence_study(argv[1]);
    return 0;
}
