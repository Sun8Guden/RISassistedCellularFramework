// Coverage probability of a user 100 m from its base station, with and without RIS-reflected signals.
#include <cstdio>
#include <risgeo/risgeo.hpp>

int main(){
    using namespace risgeo;
    const double distance = 100.0, BS_density = 1e-5, noise = 1e-13;  // m, BSs per m^2, W

    // Interfering BSs: Poisson point process beyond the serving distance, Rayleigh fading, pathloss exponent 4.
    PoissonPP2D_Reduced interference([](std::complex<double> s, double r){
        return Exponential::eval(s, PathLoss::direct_nlos(r));
    }, BS_density, distance);

    // RISs: on average 5 in the ring [10 m, 25 m] around the user, 200 elements each (CLT-Gaussian amplitude).
    const double N = 200.0;
    PoissonPP2D reflected([&](std::complex<double> s, double y, double theta){
        double variance = Parameter::antenna_gain * N * Parameter::zeta_variance
                          * PathLoss::direct_los(y) * PathLoss::reflected_los(distance, y, theta);
        double non_centrality = N * Parameter::zeta_mean * Parameter::zeta_mean / Parameter::zeta_variance;
        return ChiSquare::eval_non_central(-s, non_centrality, variance);
    }, 5.0 / (M_PI * (25.0 * 25.0 - 10.0 * 10.0)), CUBE::make_cube_2D<double>(10.0, 25.0, 0.0, 2.0 * M_PI));

    std::printf("Coverage at 0 dB: %.4f without RISs, %.4f with RISs\n",
                CoverageProbability::calculate_rayleigh(distance, 1.0, interference, noise),
                CoverageProbability::calculate_rayleigh_positive(distance, 1.0, noise, interference, reflected));
}
