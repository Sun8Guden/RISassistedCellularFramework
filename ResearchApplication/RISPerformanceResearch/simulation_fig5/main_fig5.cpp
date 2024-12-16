#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include "ReflectedSignalSquare.hpp"
#include "../StochasticGeometry/Param.hpp"
#include "../StochasticGeometry/SeparatePositive.hpp"
#include "../StochasticGeometry/ErgodicRate.hpp"
#include "../LaplaceTransformInstance/Interference.hpp"
#include "../LaplaceTransformInstance/Noise.hpp"
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../library/GenzMalik/Cube.hpp"
#include "../library/GenzMalik/GM2D.hpp"

#include "../library/json.hpp"
using json = nlohmann::json;

int main(){
    // Parameters
    double den_BS = 1e-5;
    double radius_min = 10.0; // Useless
    double radius_max = 25.0; // Useless


    double square_width = 10.0;
    double square_length = 10.0; 
    double number_RIS_panel = 3.0;
    double den_RIS = number_RIS_panel / (square_width * square_length);
    double RISnumber = 3000;
    double UEs_cell = 5.0;
    const Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);

    // Configuration 
    double distance = 100.0;
    double threshold = 1.0;

    //Experiement
    auto start = std::chrono::steady_clock::now();

    // double normative_noise_power {1e-13};

    // Define Laplace transform process for the three type RIS clusters.
    PoissonPP bs_PP(den_BS);
    // PoissonPP beamformedInterferencePP( param.beam_probability*den_RIS, cluster_cube);
    // PoissonPP scatteredInterferencePP((1.0-param.beam_probability)*den_RIS, cluster_cube);

    std::cout << " Simulation begins: " << std::endl;
    for (int i=1; i<10; i++){
        double x_ref_point = double(i) * 10.0 - square_width / 2.0;
        double y_ref_point = 0.0 - square_length / 2.0;
        auto cluster_cube = CUBE::make_cube_2D<double> (x_ref_point, (x_ref_point+square_width), y_ref_point, (y_ref_point+square_length));
        PoissonPP signalPP(den_RIS, cluster_cube);

        // double noise_power = normative_noise_power / std::pow(10.0, double(i) * 3.0 / 10.0); // varying between 0-30dB
        auto lt_interference_signal = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){
            auto lt_signal = LaplaceReflectedSignalsSquarePPP( - 1.0 * s_argument, distance_arg, cluster_cube, param);

            // Backup for BPP
            // auto lt_signal = LaplaceReflectedSignalsSquareBPP( -1.0 * s_argument, distance_arg, cluster_cube, param, number_RIS_panel );
            // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
            auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
            // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
            return  lt_signal * lt_interference; //lt_signal * lt_noise
        };
        double ergodic_rate = ergodic_rate_fix_distance(lt_interference_signal, std::complex<double>(1.0, 0.0), distance);

        std::cout << "Ergodic rate is " << ergodic_rate << std::endl;

    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() \
        << " seconds." <<std::endl;

    return 0;
}