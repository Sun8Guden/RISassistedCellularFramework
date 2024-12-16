#include <iostream>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES
#include <complex>
#include <fstream>
#include "StochasticGeometryArchiv/Param.hpp"
#include "StochasticGeometryArchiv/Cov_Prob.hpp"
#include "StochasticGeometryArchiv/Ergo_Rate.hpp"
#include "StochasticGeometryArchiv/Ergo_Rate_no_RIS.hpp"
#include "StochasticGeometryArchiv/Ergo_Rate_derivative_BS.hpp"
#include "StochasticGeometryArchiv/Ergo_Rate_derivative_RIS.hpp"
#include "library/GenzMalik/Cube.hpp"
#include "LaplaceTransformArchiv/LT_upsilon.hpp"
#include "LaplaceTransformArchiv/LT_upsilon_no_RIS.hpp"
#include "LaplaceTransformArchiv/LT_upsilon_derivative_BS.hpp"
#include "LaplaceTransformArchiv/LT_upsilon_derivative_RIS.hpp"


#include "library/json.hpp"
using json = nlohmann::json;



/* This script is to evaluate the impact of the density of BS, 

*/

int main(){
    std::string setting_direction = "./simulation_setting/";
    std::string result_direction = "./simulation_results/";
    std::string simulation_name = "parameter.json";
    
    auto start = std::chrono::steady_clock::now();

    std::ifstream ifs( (setting_direction + simulation_name) );
    json simulation_setup = json::parse(ifs);

    int num_run {1};
    std::vector<double> ergo_rate(num_run, 0.0);
    // std::vector<double> impove_BS(num_run, 0.0);
    // std::vector<double> impove_RIS(num_run, 0.0);

    for (int i = 0; i < num_run; i++){
        double den_BS = simulation_setup["den_BS"];
        den_BS *= 1e-6;
        double den_RIS = simulation_setup["den_RIS"];
        double RISnumber = simulation_setup["num_RIS_elements"];
        double radius_min = simulation_setup["radius_min"];
        double radius_max = simulation_setup["radius_max"];
        double UEs_cell = simulation_setup["UEs_cell"];
        den_RIS = den_RIS/(M_PI * (radius_max*radius_max-radius_min*radius_min));
        const Param param(den_BS, den_RIS, RISnumber, radius_min, radius_max, UEs_cell);
        auto cluster_cube = CUBE::make_cube_2D<double> (param.radius_min, param.radius_max, 0.0, 2.0 * M_PI);

        std::cout << "************PARAMETERS***************" << std::endl;
        std::cout << "Density of BS: " << simulation_setup["den_BS"] << "/km^2" << std::endl;
        std::cout << "Density of RIS: " << simulation_setup["den_RIS"] << "/cluster" << std::endl;
        std::cout << "Density of UE: " << simulation_setup["UEs_cell"] << "/cell" << std::endl;
        std::cout << "Number of RIS elements: " << simulation_setup["num_RIS_elements"] << "/RIS" << std::endl;
        std::cout << "Geometry of cluster: inner radius " << simulation_setup["radius_min"] << "m, outer radius " << simulation_setup["radius_max"] << "m" << std::endl;


        // ************** Ergodic rate with RIS. **************
        LT_upsilon lt_trans(param, cluster_cube);
        Ergo_Rate_BF ER_BR(param);
        ergo_rate.at(i) = ER_BR.calc(lt_trans, 1.0);
        
        std::cout << "************RESULTS******************" << std::endl;
        std::cout << "The ergodic rate is " <<  ergo_rate.at(i) << "nat/Hz/s" <<  std::endl;

    }

    std::cout << "************COMPUTATION TIME*********" << std::endl;
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() \
        << " seconds." <<std::endl;

    std::cout << "*************************************" << std::endl;

    json outputBS(ergo_rate);
    std::ofstream out_file( (result_direction+simulation_name) );
    out_file << outputBS << std::endl;

    return 0; 
}