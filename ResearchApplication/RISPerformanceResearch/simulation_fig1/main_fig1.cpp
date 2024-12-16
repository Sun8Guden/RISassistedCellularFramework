#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <complex>
#include <fstream>
#include "../StochasticGeometry/Param.hpp"
#include "../StochasticGeometry/SeparatePositive.hpp"
#include "../StochasticGeometry/ErgodicRate.hpp"
#include "../LaplaceTransformInstance/Interference.hpp"
#include "../LaplaceTransformInstance/Noise.hpp"
#include "../LaplaceTransformInstance/ReflectedSignal.hpp"
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../library/GenzMalik/Cube.hpp"
#include "../library/GenzMalik/GM2D.hpp"

#include "../library/json.hpp"
using json = nlohmann::json;

int main(){
    // Parameters
    double radius_min = 10.0;
    double radius_max = 25.0;
    double den_RIS = 5.0 / (M_PI * (radius_max*radius_max-radius_min*radius_min));
    double RISnumber = 3000;
    double UEs_cell = 5.0;
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    // Configuration 
    double distance = 100.0;
    double threshold = 1.0;

    //Experiement
    auto start = std::chrono::steady_clock::now();
    double normative_noise_power {1e-13};

    // Define Laplace transform process for the three type RIS clusters.
    PoissonPP signalPP(den_RIS, cluster_cube);

    
    double noise_power = normative_noise_power / std::pow(10.0, 0.0); // 0dB-30dB, 

    std::vector<double> rate_SINR_RIS;
    std::vector<double> rate_SIR_RIS;
    std::vector<double> rate_SINR_BS;
    std::vector<double> rate_SIR_BS;
    std::vector<double> rate_SNR;

    for (int i=1; i<=10; i++){
        double den_BS = 1e-6 * double(i);
        PoissonPP bs_PP(den_BS);
        const Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);
        PoissonPP beamformedInterferencePP( param.beam_probability*den_RIS, cluster_cube);
        PoissonPP scatteredInterferencePP((1.0-param.beam_probability)*den_RIS, cluster_cube);

    
        auto lt_SINR_RIS = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

            auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
            auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
            // auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
            auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
            return  lt_signal * lt_interference * lt_noise; //lt_signal 
        };
        double ergodic_rate_SINR_RIS = ergodic_rate_fix_distance(lt_SINR_RIS, std::complex<double>(1.0, 0.0), distance);
        rate_SINR_RIS.push_back(ergodic_rate_SINR_RIS);

        // auto lt_SIR_RIS = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

        //     auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
        //     auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
        //     // auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
        //     // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
        //     return  lt_signal * lt_interference;// * lt_noise; //lt_signal 
        // };
        // double ergodic_rate_SIR_RIS = ergodic_rate_fix_distance(lt_SIR_RIS, std::complex<double>(1.0, 0.0), distance);
        // rate_SIR_RIS.push_back(ergodic_rate_SIR_RIS);

        // auto lt_SINR_BS = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

        //     auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
        //     // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
        //     auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
        //     auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
        //     return  lt_signal * lt_interference * lt_noise; //lt_signal 
        // };
        // double ergodic_rate_SINR_BS = ergodic_rate_fix_distance(lt_SINR_BS, std::complex<double>(1.0, 0.0), distance);
        // rate_SINR_BS.push_back(ergodic_rate_SINR_BS);

        // auto lt_SIR_BS = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

        //     auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
        //     // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
        //     auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
        //     // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
        //     return  lt_signal * lt_interference;// * lt_noise; //lt_signal 
        // };
        // double ergodic_rate_SIR_BS = ergodic_rate_fix_distance(lt_SIR_BS, std::complex<double>(1.0, 0.0), distance);
        // rate_SIR_BS.push_back(ergodic_rate_SIR_BS);        

        // auto lt_SNR = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg){

        //     auto lt_signal = LaplaceReflectedSignals( - 1.0 * s_argument, distance_arg, param, signalPP);
        //     // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
        //     // auto lt_interference = LT_interference_BSs(s_argument, distance_arg, threshold_arg, param, bs_PP);  // Interference from BS
        //     auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
        //     return  lt_signal * lt_noise; //lt_signal 
        // };
        // double ergodic_rate_SNR = ergodic_rate_fix_distance(lt_SNR, std::complex<double>(1.0, 0.0), distance);
        // rate_SNR.push_back(ergodic_rate_SNR);

    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() \
        << " seconds." <<std::endl;


    json output_rate_SINR_RIS( rate_SINR_RIS );
    // json output_rate_SIR_RIS( rate_SIR_RIS );
    // json output_rate_SINR_BS( rate_SINR_BS );
    // json output_rate_SIR_BS( rate_SIR_BS );
    // json output_rate_SNR( rate_SNR );

    std::ofstream out_file( "output0dBnoise.json" );

    
    out_file << "output_rate_SINR_RIS = " << output_rate_SINR_RIS << ";" <<std::endl;
    // out_file << "output_rate_SIR_RIS = " << output_rate_SIR_RIS << ";" <<std::endl;
    // out_file << "output_rate_SINR_BS = " << output_rate_SINR_BS << ";" <<std::endl;
    // out_file << "output_rate_SIR_BS = " << output_rate_SIR_BS << ";" <<std::endl;
    // out_file << "output_rate_SNR = " << output_rate_SNR << ";" <<std::endl;

    return 0;
}