#include <iostream>
#include <fstream>
#include <cmath>
#include "Interference_mean.hpp"

#include "../../../Utils/json.hpp"
using json = nlohmann::json;



int main(){
    double den_BS = 1e-5;
    double radius_min = 10.0;
    double radius_max = 25.0;
    double den_RIS = 5.0 / (M_PI * (radius_max*radius_max-radius_min*radius_min));
    double RISnumber = 3000;
    double UEs_cell = 5.0;
    Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);
    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);

    // Configuration 
    double distance = 100.0;
    double threshold = 1.0;

    //Experiement
    auto start = std::chrono::steady_clock::now();
    PoissonPP signalPP(den_RIS, cluster_cube);
    PoissonPP bs_PP(den_BS);

    std::vector<double> Interference_BS;
    std::vector<double> Interference_RIS_beam;
    std::vector<double> Interference_RIS_scatter;
    std::vector<double> Interference_RIS;
    std::vector<double> Interference_all;


    for (int i=0; i <= 36; i=i+2){
        // param.beam_probability = 0.5 / (std::pow( 10.0, 2.0 * ( double(i) / 10.0 )));
        param.beam_probability =  double(i) / 360.0; // few degrees out of 360 degrees
        PoissonPP beamformedInterferencePP( param.beam_probability * den_RIS, cluster_cube);
        PoissonPP scatteredInterferencePP( ( 1.0 - param.beam_probability ) * den_RIS, cluster_cube);
  
        // Metrics
        double cur_average_interference_BS = LT_interference_mean_BS(0.0, distance, threshold, param, bs_PP);
        Interference_BS.push_back(cur_average_interference_BS);
        double cur_average_interference_RIS_beam = LT_interference_mean_beamform(0.0, distance, threshold, param, bs_PP, beamformedInterferencePP);
        Interference_RIS_beam.push_back(cur_average_interference_RIS_beam);
        double cur_average_interference_RIS_scatter = LT_interference_mean_scatter(0.0, distance, threshold, param, bs_PP, scatteredInterferencePP);
        Interference_RIS_scatter.push_back(cur_average_interference_RIS_scatter);
        double cur_average_interference_RIS = LT_interference_mean_RIS(0.0, distance, threshold, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP);
        Interference_RIS.push_back(cur_average_interference_RIS);
        // std::cout << "*************************** RIS discussion ***************************" << std::endl;
        // std::cout << "Beam part is " << cur_average_interference_RIS_beam / cur_average_interference_RIS << std::endl;
        // std::cout << "Scattered part is " << cur_average_interference_RIS_scatter / cur_average_interference_RIS << std::endl;
        double cur_average_interference_all = LT_interference_mean_BS_RIS(0.0, distance, threshold, param, bs_PP, beamformedInterferencePP, scatteredInterferencePP);
        Interference_all.push_back(cur_average_interference_all);
        // std::cout << "*************************** Interference discussion ***************************" << std::endl;
        // std::cout << "RIS part is " << cur_average_interference_RIS/cur_average_interference_all << std::endl;
        // std::cout << "BS part is " << cur_average_interference_BS/cur_average_interference_all << std::endl;
        double ratio_of_RIS = cur_average_interference_RIS/cur_average_interference_all;
    }


    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() \
        << " microseconds." <<std::endl;

    json output_BS( Interference_BS );
    json output_RIS( Interference_RIS );
    json output_RIS_beam( Interference_RIS_beam );
    json output_RIS_scatter( Interference_RIS_scatter );
    json output_ALL( Interference_all );

    std::ofstream out_file( "output.json" );
    out_file << "Interference_BS=" << output_BS <<";" << std::endl;
    out_file << "Interference_RIS=" << output_RIS <<";" << std::endl;
    out_file << "Interference_RIS_beam=" << output_RIS_beam <<";" << std::endl;
    out_file << "Interference_RIS_scatter=" << output_RIS_scatter <<";" << std::endl;
    out_file << "Interference_all=" << output_ALL <<";" << std::endl;



    // std::cout << "Power of interferece from BSs: " << average_interference_BS << std::endl;
    // std::cout << "Power of interferece from RISs: " << average_interference_RIS << std::endl;
    // std::cout << "Power of interferece from BSs and RISs: " << average_interference_all << std::endl;

}