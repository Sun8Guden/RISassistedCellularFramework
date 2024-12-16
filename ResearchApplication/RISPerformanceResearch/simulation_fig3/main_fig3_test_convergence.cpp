#include <iostream>
#include <cmath>
#include <chrono>
#include <complex>
#include <fstream>
#include "../StochasticGeometry/Param.hpp"
#include "../StochasticGeometry/SeparatePositive.hpp"
#include "../StochasticGeometry/ErgodicRate.hpp"
#include "ErgodicRateCoverageHole.hpp"
#include "../LaplaceTransformInstance/Interference.hpp"
#include "../LaplaceTransformInstance/Noise.hpp"
// #include "../LaplaceTransformInstance/ReflectedSignal.hpp"
#include "ReflectedSignalCoverageHole.hpp"
#include "../LaplaceTransformPointProcess/PoissonPP.hpp"
#include "../library/GenzMalik/Cube.hpp"
#include "../library/GenzMalik/GM2D.hpp"

#include "../library/json.hpp"
using json = nlohmann::json;

int main(){
    // Parameters
    double den_BS = 4e-6;
    double radius_min = 25.0;
    double radius_max = 30.0;
    double radius_hole = 10.0;

    double den_RIS = 4.0 / (M_PI * (radius_max*radius_max-radius_min*radius_min));
    double RISnumber = 3000;
    double UEs_cell = 5.0;
    const Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    auto cluster_UE_ball = CUBE::make_cube_2D<double> (0.0, radius_hole, 0.0, 2.0*M_PI);

    // Configuration 
    double distance = 80.0;
    double threshold = 1.0;

    //Experiement
    auto start = std::chrono::steady_clock::now();

    double normative_noise_power {1e-13};

    // Define Laplace transform process for the three type RIS clusters.
    PoissonPP signalPP(den_RIS, cluster_cube);
    PoissonPP bs_PP(den_BS);
    // PoissonPP beamformedInterferencePP( param.beam_probability*den_RIS, cluster_cube);
    // PoissonPP scatteredInterferencePP((1.0-param.beam_probability)*den_RIS, cluster_cube);


    std::vector<double> ergo_rate_log;
    for (int i=0; i<=10; i++){
        double panelty = std::pow(10.0, double(i) / 10.0);

        // std::cout << "The panelty is " << panelty  <<std::endl;
        // auto lt_interference_signal = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg, double UE_hole_radius, double UE_hole_angle){
        //     auto lt_signal = LaplaceReflectedSignalsCoverageHolePPP(- 1.0 * s_argument, distance_arg, panelty, cluster_cube, param, UE_hole_radius, UE_hole_angle);
        //     // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
        //     auto lt_interference = LT_interference_BSs(s_argument, distance_arg, panelty * threshold_arg, param, bs_PP);  // Interference from BS
        //     // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
        //     return  lt_signal * lt_interference; // lt_noise
        // };
        // double ergodic_rate = ergodic_rate_coverage_hole(lt_interference_signal, std::complex<double>(1.0, 0.0), distance, radius_hole);


        auto lt_interference_signal = [&](std::complex<double> s_argument, double distance_arg, double threshold_arg, double x_UE, double y_UE){

            auto lt_signal = LaplaceReflectedSignalsCoverageHoleFix( - 1.0 * s_argument, distance_arg, panelty, cluster_cube, param, x_UE, y_UE);
            // auto lt_interference = LT_interference(s_argument, distance_arg, threshold_arg, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP); // Interference from BS + RISs
            auto lt_interference = LT_interference_BSs(s_argument, distance_arg, panelty * threshold_arg, param, bs_PP);  // Interference from BS
            // auto lt_noise = LT_noise(s_argument, noise_power, distance_arg, threshold_arg, param);
            return  lt_signal * lt_interference; //   lt_noise
        };
        

        auto ergodic_rate_area = [&](double r_UE, double theta_UE){
            double x_UE = r_UE * std::cos(theta_UE);
            double y_UE = r_UE * std::sin(theta_UE);
            double cur_ergodic_rate = ergodic_rate_fix_distance_coverage_hole(lt_interference_signal, std::complex<double>(1.0, 0.0), distance, x_UE, y_UE);
            return cur_ergodic_rate;
        };

        // double ergodic_rate = ergodic_rate_area(10.0, 0.0);

        double estimated_error{0.0};    
        unsigned int number_evaluation{0};
        double ergodic_rate = ergodic_rate_area(10.0, 0);
        // double ergodic_rate = GM::GM2D<double>::integrate(ergodic_rate_area, cluster_UE_ball, 1e-3, estimated_error, 1000, number_evaluation);
        // ergodic_rate = ergodic_rate / (M_PI * radius_hole*radius_hole);
        std::cout << "Ergodic rate is " << ergodic_rate << ", with function call " << \
            number_evaluation << " and estimated error " << estimated_error <<   std::endl;
        ergo_rate_log.push_back(ergodic_rate);
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() \
        << " seconds." <<std::endl;


    json output_rate_log( ergo_rate_log );
    std::ofstream out_file( "output.json" );
    out_file << output_rate_log << std::endl;

    return 0;
}