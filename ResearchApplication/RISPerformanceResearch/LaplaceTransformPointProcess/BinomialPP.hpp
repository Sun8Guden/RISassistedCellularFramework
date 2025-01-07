#pragma once
#include <iostream>
#include <cassert>
#include <complex>
#include <utility>
#include <cmath>
#include "../../../Integration/GenzMalik/Cube.hpp"
#include "../../../Integration/GenzMalik/GM2D.hpp"

class BinomialPP {

    int number_points;
    CUBE::Cube<double, 2> cluster_cube;
    double probability;


public:
    BinomialPP(int number_points_, const CUBE::Cube<double, 2>& cluster_cube_): number_points(number_points_), cluster_cube(cluster_cube_) {
        double radius_max = cluster_cube.cube_center.at(0) + 0.5 * cluster_cube.cube_size.at(0);
        double radius_min = cluster_cube.cube_center.at(0) - 0.5 * cluster_cube.cube_size.at(0);
        probability = 1.0 / ( M_PI * (radius_max * radius_max - radius_min * radius_min) );
    }
    BinomialPP(int number_points_) : number_points(number_points_){}
    ~BinomialPP(){}


    template<class F>
    auto eval_2D(F f)->decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) const
    {
        typedef decltype(std::declval<F>()(std::declval<double>(), std::declval<double>())) return_type;
        double estimated_error{0.0};    
        unsigned int number_evaluation{0};

        return_type lt_value = GM::GM2D<double>::integrate(f, cluster_cube, 1e-6, estimated_error, 5000, number_evaluation);
        // std::cout << probability * lt_value << "probability "<< probability <<"lt_value " <<lt_value << std::endl; 
        return std::pow( probability * lt_value,  double( number_points ) );
    }

};