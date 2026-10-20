// Monte Carlo validation of the analytical coverage probability and ergodic rate.
//
// Simulates the network of the examples (and of validation/matlab): interfering BSs form a PPP outside the
// serving distance, all BS links are Rayleigh faded with NLoS pathloss, and a Poisson number of RISs is
// uniformly placed in a ring around the user, each reflecting a CLT-Gaussian amplitude over two LoS hops.
// The analytical results from the library are compared with the simulation at every SINR threshold.
//
// Samples are drawn in fixed chunks with their own seeds, so results do not depend on the number of threads.
//
// Usage: monte_carlo [samples] [output_dir]     (defaults: 1000000 samples, no CSV output)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "../include/risgeo/risgeo.hpp"

using namespace risgeo;

namespace {

// Scenario of examples/ and validation/matlab/.
constexpr double noise_power = 1e-13;          // W
constexpr double BS_density = 1e-5;            // per m^2
constexpr double distance_typical_BS = 100.0;  // m
constexpr double simulation_radius = 3000.0;   // m, interferers beyond are negligible
constexpr double mean_RIS_per_cluster = 5.0;
constexpr double RIS_elements = 200.0;
constexpr double radius_min = 10.0;            // m, RIS ring around the user
constexpr double radius_max = 25.0;

double seconds_since(std::chrono::steady_clock::time_point start){
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

struct SinrSample {
    double without_RIS;
    double with_RIS;
};

class NetworkSampler {
    std::mt19937_64 m_rng;
    std::poisson_distribution<int> m_number_BS;
    std::poisson_distribution<int> m_number_RIS;
    std::exponential_distribution<double> m_rayleigh_power{1.0};
    std::uniform_real_distribution<double> m_uniform{0.0, 1.0};
    std::normal_distribution<double> m_gaussian{0.0, 1.0};

    // Uniform point in the ring [r_min, r_max] (in the plane), returns its radius.
    double ring_radius(double r_min, double r_max){
        return std::sqrt(m_uniform(m_rng) * (r_max * r_max - r_min * r_min) + r_min * r_min);
    }

public:
    explicit NetworkSampler(std::uint64_t seed)
        : m_rng(seed),
          m_number_BS(BS_density * M_PI * (simulation_radius * simulation_radius - distance_typical_BS * distance_typical_BS)),
          m_number_RIS(mean_RIS_per_cluster) {}

    SinrSample draw(){
        const double gain = Parameter::antenna_gain * Parameter::transmission_power;

        double interference = 0.0;
        const int number_BS = m_number_BS(m_rng);
        for (int i = 0; i < number_BS; i++) {
            interference += m_rayleigh_power(m_rng) * gain * PathLoss::direct_nlos(ring_radius(distance_typical_BS, simulation_radius));
        }
        const double signal_BS = m_rayleigh_power(m_rng) * gain * PathLoss::direct_nlos(distance_typical_BS);

        // Each RIS reflects the sum of RIS_elements random element gains: Gaussian amplitude by the CLT.
        const double amplitude_mean = Parameter::zeta_mean * RIS_elements;
        const double amplitude_std = std::sqrt(Parameter::zeta_variance * RIS_elements);
        double signal_RIS = 0.0;
        const int number_RIS = m_number_RIS(m_rng);
        for (int j = 0; j < number_RIS; j++) {
            const double y = ring_radius(radius_min, radius_max);
            const double theta = 2.0 * M_PI * m_uniform(m_rng);
            const double amplitude = amplitude_mean + amplitude_std * m_gaussian(m_rng);
            // BS -> RIS and RIS -> user are both LoS links with their own antenna gain.
            signal_RIS += amplitude * amplitude * gain * PathLoss::reflected_los(distance_typical_BS, y, theta)
                          * Parameter::antenna_gain * PathLoss::direct_los(y);
        }

        const double denominator = interference + noise_power;
        return {signal_BS / denominator, (signal_BS + signal_RIS) / denominator};
    }
};

struct Tally {
    std::vector<double> covered_without_RIS, covered_with_RIS;
    double rate_without_RIS{0.0}, rate_with_RIS{0.0};
    double rate_without_RIS_sq{0.0}, rate_with_RIS_sq{0.0};

    explicit Tally(std::size_t thresholds) : covered_without_RIS(thresholds, 0.0), covered_with_RIS(thresholds, 0.0) {}

    void add(const Tally& other){
        for (std::size_t k = 0; k < covered_without_RIS.size(); k++) {
            covered_without_RIS[k] += other.covered_without_RIS[k];
            covered_with_RIS[k] += other.covered_with_RIS[k];
        }
        rate_without_RIS += other.rate_without_RIS;
        rate_with_RIS += other.rate_with_RIS;
        rate_without_RIS_sq += other.rate_without_RIS_sq;
        rate_with_RIS_sq += other.rate_with_RIS_sq;
    }
};

double standard_error_proportion(double p, double n){ return std::sqrt(p * (1.0 - p) / n); }

// z-score of the simulated proportion under the hypothesis that the analytical value is exact.
double z_score(double analytic, double simulated, double n){
    return (simulated - analytic) / std::max(standard_error_proportion(analytic, n), 1.0 / n);
}

} // namespace

int main(int argc, char** argv){
    const long samples = argc > 1 ? std::atol(argv[1]) : 1000000;
    const std::string output_dir = argc > 2 ? argv[2] : "";
    const long chunk = 10000;
    const long chunks = (samples + chunk - 1) / chunk;
    const std::uint64_t seed = 20231010;

    std::vector<double> thresholds_dB;
    for (double t = -10.0; t <= 20.0 + 1e-9; t += 2.5) thresholds_dB.push_back(t);
    std::vector<double> thresholds;
    for (double t : thresholds_dB) thresholds.push_back(std::pow(10.0, t / 10.0));

    // ---------------- Monte Carlo ----------------
    auto start_mc = std::chrono::steady_clock::now();
    std::vector<Tally> chunk_tallies(chunks, Tally(thresholds.size()));
    #pragma omp parallel for schedule(dynamic)
    for (long c = 0; c < chunks; c++) {
        NetworkSampler sampler(seed + std::uint64_t(c));
        Tally& tally = chunk_tallies[c];
        const long n = std::min(chunk, samples - c * chunk);
        for (long i = 0; i < n; i++) {
            const SinrSample sinr = sampler.draw();
            for (std::size_t k = 0; k < thresholds.size(); k++) {
                tally.covered_without_RIS[k] += sinr.without_RIS > thresholds[k];
                tally.covered_with_RIS[k] += sinr.with_RIS > thresholds[k];
            }
            const double r0 = std::log1p(sinr.without_RIS), r1 = std::log1p(sinr.with_RIS);
            tally.rate_without_RIS += r0;  tally.rate_without_RIS_sq += r0 * r0;
            tally.rate_with_RIS += r1;     tally.rate_with_RIS_sq += r1 * r1;
        }
    }
    Tally total(thresholds.size());
    for (const Tally& t : chunk_tallies) total.add(t);
    const double mc_seconds = seconds_since(start_mc);
    const double n = double(samples);

    // ---------------- Analytical ----------------
    auto cluster_cube = CUBE::make_cube_2D<double>(radius_min, radius_max, 0.0, 2.0 * M_PI);
    const double RIS_density = mean_RIS_per_cluster / (M_PI * (radius_max * radius_max - radius_min * radius_min));
    // The interferers are integrated over the same region as simulated, [distance, simulation_radius];
    // integrating to infinity instead lowers the coverage by about 0.2 % (relative) at 7.5 dB.
    PoissonPP2D_Reduced interference([](std::complex<double> s, double r){
        return Exponential::eval(s, PathLoss::direct_nlos(r));
    }, BS_density, distance_typical_BS, simulation_radius);
    PoissonPP2D reflected_signals([](std::complex<double> s, double y, double theta){
        double variance = Parameter::antenna_gain * RIS_elements * Parameter::zeta_variance
                          * PathLoss::direct_los(y) * PathLoss::reflected_los(distance_typical_BS, y, theta);
        double non_centrality = RIS_elements * Parameter::zeta_mean * Parameter::zeta_mean / Parameter::zeta_variance;
        return ChiSquare::eval_non_central(-s, non_centrality, variance);
    }, RIS_density, cluster_cube);

    std::vector<double> analytic_without_RIS, analytic_with_RIS;
    auto start_analytic = std::chrono::steady_clock::now();
    for (double t : thresholds) {
        analytic_without_RIS.push_back(CoverageProbability::calculate_rayleigh(distance_typical_BS, t, interference, noise_power));
    }
    const double analytic_without_seconds = seconds_since(start_analytic) / thresholds.size();
    start_analytic = std::chrono::steady_clock::now();
    for (double t : thresholds) {
        analytic_with_RIS.push_back(CoverageProbability::calculate_rayleigh_positive(distance_typical_BS, t, noise_power, interference, reflected_signals));
    }
    const double analytic_with_seconds = seconds_since(start_analytic) / thresholds.size();
    start_analytic = std::chrono::steady_clock::now();
    const double rate_without_RIS = ErgodicRate::calculate_rayleigh(distance_typical_BS, interference, noise_power);
    const double rate_without_seconds = seconds_since(start_analytic);
    start_analytic = std::chrono::steady_clock::now();
    const double rate_with_RIS = ErgodicRate::calculate_rayleigh_positive(distance_typical_BS, interference, noise_power, reflected_signals);
    const double rate_with_seconds = seconds_since(start_analytic);

    // ---------------- Report ----------------
    std::printf("Monte Carlo: %ld samples in %.2f s (%.0f ns/sample)\n\n", samples, mc_seconds, 1e9 * mc_seconds / n);
    std::printf("%10s | %10s %10s %7s | %10s %10s %7s\n", "SINR [dB]", "analytic", "MC", "z", "analytic", "MC", "z");
    std::printf("%10s | %29s | %29s\n", "", "coverage without RIS", "coverage with RIS");
    double max_z = 0.0;
    std::vector<double> mc_without(thresholds.size()), mc_with(thresholds.size());
    for (std::size_t k = 0; k < thresholds.size(); k++) {
        mc_without[k] = total.covered_without_RIS[k] / n;
        mc_with[k] = total.covered_with_RIS[k] / n;
        double z0 = z_score(analytic_without_RIS[k], mc_without[k], n);
        double z1 = z_score(analytic_with_RIS[k], mc_with[k], n);
        max_z = std::max({max_z, std::abs(z0), std::abs(z1)});
        std::printf("%10.1f | %10.6f %10.6f %7.2f | %10.6f %10.6f %7.2f\n", thresholds_dB[k],
                    analytic_without_RIS[k], mc_without[k], z0, analytic_with_RIS[k], mc_with[k], z1);
    }
    auto rate_se = [&](double sum, double sum_sq){ double m = sum / n; return std::sqrt((sum_sq / n - m * m) / n); };
    const double mc_rate_without = total.rate_without_RIS / n, mc_rate_with = total.rate_with_RIS / n;
    const double se_without = rate_se(total.rate_without_RIS, total.rate_without_RIS_sq);
    const double se_with = rate_se(total.rate_with_RIS, total.rate_with_RIS_sq);
    std::printf("\nErgodic rate [nat/s/Hz]   analytic   MC (+- 1 s.e.)\n");
    std::printf("  without RIS            %9.5f   %9.5f +- %.5f\n", rate_without_RIS, mc_rate_without, se_without);
    std::printf("  with RIS               %9.5f   %9.5f +- %.5f\n", rate_with_RIS, mc_rate_with, se_with);
    std::printf("\nLargest |MC - analytic| / standard error over all thresholds: %.2f\n", max_z);
    std::printf("Analytic time per point: coverage %.2e s, coverage with RIS %.2e s, rate %.2e s, rate with RIS %.2f s\n",
                analytic_without_seconds, analytic_with_seconds, rate_without_seconds, rate_with_seconds);

    if (!output_dir.empty()) {
        std::ofstream coverage(output_dir + "/coverage_validation.csv");
        coverage << "threshold_db,analytic_without_ris,mc_without_ris,mc_without_ris_se,analytic_with_ris,mc_with_ris,mc_with_ris_se\n";
        for (std::size_t k = 0; k < thresholds.size(); k++) {
            coverage << thresholds_dB[k] << "," << analytic_without_RIS[k] << "," << mc_without[k] << ","
                     << standard_error_proportion(mc_without[k], n) << "," << analytic_with_RIS[k] << "," << mc_with[k] << ","
                     << standard_error_proportion(mc_with[k], n) << "\n";
        }

        // Accuracy vs. cost at 0 dB with RISs, single thread: observed error of the running Monte Carlo estimate
        // against the analytical value, and the cost of the analytical evaluation itself.
        std::ofstream convergence(output_dir + "/accuracy_vs_time.csv");
        convergence << "method,samples,seconds,absolute_error,probability\n";
        const std::size_t k0 = 4;  // 0 dB
        const double reference = analytic_with_RIS[k0];
        NetworkSampler sampler(seed + 999999);
        auto start = std::chrono::steady_clock::now();
        double covered = 0.0;
        long next_checkpoint = 100;
        for (long i = 1; i <= samples; i++) {
            covered += sampler.draw().with_RIS > thresholds[k0];
            if (i == next_checkpoint || i == samples) {
                convergence << "monte_carlo," << i << "," << seconds_since(start) << "," << std::abs(covered / i - reference)
                            << "," << reference << "\n";
                next_checkpoint = long(next_checkpoint * 1.25) + 1;
            }
        }
        convergence << "analytic,," << analytic_with_seconds << "," << Parameter::expected_error << "," << reference << "\n";

        std::ofstream summary(output_dir + "/validation_summary.csv");
        summary << "quantity,analytic,monte_carlo,monte_carlo_se,analytic_seconds\n";
        summary << "ergodic_rate_without_ris," << rate_without_RIS << "," << mc_rate_without << "," << se_without << "," << rate_without_seconds << "\n";
        summary << "ergodic_rate_with_ris," << rate_with_RIS << "," << mc_rate_with << "," << se_with << "," << rate_with_seconds << "\n";
        summary << "monte_carlo_seconds_per_sample,,," << "," << mc_seconds / n << "\n";
        std::printf("CSV files written to %s\n", output_dir.c_str());
    }
    return max_z < 5.0 ? 0 : 1;
}
