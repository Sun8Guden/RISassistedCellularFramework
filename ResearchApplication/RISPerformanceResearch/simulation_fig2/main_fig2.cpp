#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include <fstream>
#include "../StochasticGeometryHelper/Param.hpp"
#include "../StochasticGeometryHelper/SeparatePositive.hpp"
#include "../StochasticGeometryHelper/ErgodicRate.hpp"
#include "../LaplaceTransformFunctions/Interference.hpp"
#include "../LaplaceTransformFunctions/Noise.hpp"
#include "../LaplaceTransformFunctions/ReflectedSignal.hpp"
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"

#include "../../../Utils/json.hpp"
using json = nlohmann::json;

int main(){
    // Parameters
    double den_BS = 1e-5;
    double radius_min = 15.0;
    double radius_max = 35.0;
    double UEs_cell = 5.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    // Configuration 
    double distance = 100.0;
    double threshold = 1.0;

    auto start = std::chrono::steady_clock::now();

    double normative_noise_power {1e-13};
    std::vector<double> ergo_rate_vec;

    for (int i=8; i>=2; i--){
        double den_RIS = double(i) / (M_PI * (radius_max*radius_max-radius_min*radius_min));
        double total_RIS_number = 10000.0;
        double RISnumber = total_RIS_number/double(i);
        const Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);

        // Define Laplace transform process for the three type RIS clusters.
        PoissonPP signalPP(den_RIS, cluster_cube);
        PoissonPP bs_PP(den_BS);
        // PoissonPP beamformedInterferencePP( param.beam_probability*den_RIS, cluster_cube);
        // PoissonPP scatteredInterferencePP((1.0-param.beam_probability)*den_RIS, cluster_cube);


        // double noise_power = normative_noise_power / std::pow(10.0, double(i) * 3.0 / 10.0); // varying between 0-30dB
        auto lt_interference_signal = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

            auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
            // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
            auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
            // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
            return  lt_signal * lt_interference; //lt_signal * lt_noise
        };
        double ergodic_rate = ergodic_rate_fix_distance(lt_interference_signal, std::complex<double>(1.0, 0.0), distance);
        ergo_rate_vec.push_back(ergodic_rate);
        // std::cout << "Ergodic rate is " << ergodic_rate << std::endl;

    }

    json output_rate( ergo_rate_vec );
    std::ofstream out_file( "output.json" );
    out_file << "output_rate = " << output_rate << ";" <<std::endl;


    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() \
        << " seconds." <<std::endl;

    return 0;
}