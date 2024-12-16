#include <iostream>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES
#include <complex>
#include <fstream>
#include "StochasticGeometry/Param.hpp"
#include "StochasticGeometry/Cov_Prob.hpp"
#include "StochasticGeometry/Ergo_Rate.hpp"
#include "StochasticGeometry/Ergo_Rate_no_RIS.hpp"
#include "StochasticGeometry/Ergo_Rate_derivative_BS.hpp"
#include "StochasticGeometry/Ergo_Rate_derivative_RIS.hpp"
#include "library/GenzMalik/Cube.hpp"
#include "Laplace_transform/LT_upsilon.hpp"
#include "Laplace_transform/LT_upsilon_no_RIS.hpp"
#include "Laplace_transform/LT_upsilon_derivative_BS.hpp"
#include "Laplace_transform/LT_upsilon_derivative_RIS.hpp"
#include "library/json.hpp"
using json = nlohmann::json;

int main(){

    double Cost_BS_RIS = 10.0;

    std::string setting_direction = "./simulation_setting/";
    std::string result_direction = "./simulation_results/";
    std::string simulation_name = "investment.json";
    std::ifstream ifs( (setting_direction + simulation_name) );
    json simulation_setup = json::parse(ifs);

    double den_BS = simulation_setup["den_BS"];
    double den_RIS = simulation_setup["den_RIS"];
    double RISnumber = simulation_setup["num_RIS_elements"];
    double radius_min = simulation_setup["radius_min"];
    double radius_max = simulation_setup["radius_max"];
    double UEs_cell = simulation_setup["UEs_cell"];

    double BS_den_init = den_BS * 1e-6;
    double RIS_den_init = den_RIS;

    double improvement_by_BS{0.0};
    double improvement_by_RIS{0.0};

    int investment_round {50};
    double step_size = 0.2;

    std::vector<double> capacity_vec(investment_round, 0.0);
    std::vector<double> BS_den_vec(investment_round, 0.0);
    std::vector<double> RIS_den_vec(investment_round, 0.0);
    std::vector<double> BS_gain_vec(investment_round, 0.0);
    std::vector<double> RIS_gain_vec(investment_round, 0.0);
    std::vector<bool> bet_on_BS(investment_round, true); 

    auto cluster_cube = CUBE::make_cube_2D<double> (radius_min, radius_max, 0.0, 2.0 * M_PI);
    bool UsingBS {true};

    for (int i=0; i < investment_round; i++){

        double den_RIS_cluster = RIS_den_init / BS_den_init/ (M_PI * (radius_max*radius_max-radius_min*radius_min));
        const Param param(BS_den_init, den_RIS_cluster, RISnumber, radius_min, radius_max, UEs_cell);

        LT_upsilon lt_trans(param, cluster_cube);
        LT_upsilon_derivative_BS lt_BS(param, cluster_cube);
        LT_upsilon_derivative_RIS lt_RIS(param, cluster_cube);

        Ergo_Rate_BF er_BR(param);
        Ergo_Rate_BS_derivative er_BS_improve(param);
        Ergo_Rate_RIS_derivative er_RIS_improve(param);

        if (den_RIS_cluster == 0.0){
            improvement_by_BS = er_BS_improve.calc_no_RIS(lt_BS, 1.0); //correct
            improvement_by_RIS = er_RIS_improve.calc_no_RIS(lt_RIS, 1.0);
            capacity_vec.at(i) = er_BR.calc_no_RIS(lt_trans, 1.0);
            improvement_by_RIS = improvement_by_RIS / BS_den_init / (M_PI * (radius_max*radius_max-radius_min*radius_min));
        } else {
            improvement_by_BS = er_BS_improve.calc(lt_BS, 1.0);
            improvement_by_RIS = er_RIS_improve.calc(lt_RIS, 1.0);
            improvement_by_RIS = improvement_by_RIS / BS_den_init / (M_PI * (radius_max*radius_max-radius_min*radius_min));
            capacity_vec.at(i) = er_BR.calc(lt_trans, 1.0);
        }
        if (improvement_by_BS/Cost_BS_RIS > improvement_by_RIS){
            bet_on_BS.at(i) = true;
            BS_den_init += step_size * 1e-6;
        } else {
            // UsingBS = false;
            bet_on_BS.at(i) = false;
            RIS_den_init += step_size * Cost_BS_RIS * 1e-6; // density per m^2
        }

        BS_den_vec.at(i) = BS_den_init;
        RIS_den_vec.at(i) = RIS_den_init/BS_den_init; // Number of RIS per cluster. 
        BS_gain_vec.at(i) = improvement_by_BS * 1e-6;
        RIS_gain_vec.at(i) = improvement_by_RIS * Cost_BS_RIS * 1e-6;

    }

    json output_rate(capacity_vec);
    json output_BS(BS_den_vec);
    json output_RIS(RIS_den_vec);
    json output_decision(bet_on_BS);
    json output_BS_gain(BS_gain_vec);
    json output_RIS_gain(RIS_gain_vec);

    std::ofstream out_file( (result_direction+simulation_name) );
    out_file << "Ergodic rate: " << output_rate << std::endl;
    out_file << "Density BS: " << output_BS << std::endl;
    out_file << "Density RIS: " <<  output_RIS << std::endl;
    out_file << "slope BS: " << output_BS_gain << std::endl;
    out_file << "slope RIS: " <<output_RIS_gain << std::endl;
    out_file << "Bet on BS?: "<< output_decision << std::endl;
   

    return 0;
}