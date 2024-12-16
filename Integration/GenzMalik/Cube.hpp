#pragma once

#include <array>
#include <vector>
#include <stdexcept>

namespace CUBE{

    template <typename Real, int N>
    class Cube {
    public:
        std::array<Real, N> cube_center;
        std::array<Real, N> cube_size;

        Real esti_error{ 0.0 };
        
        bool evaluated{ false };
        unsigned int split_direction{ 0 };

        Cube(){}
        ~Cube(){}

        Cube(const std::array<Real, N>& center_arg, const std::array<Real, N>& size_arg) : \
            cube_center(center_arg), cube_size(size_arg){}

        // Copy constructor
        Cube(const Cube& source_cube){
            cube_center = source_cube.cube_center;
            cube_size = source_cube.cube_size;
            evaluated = source_cube.evaluated;
            esti_error = source_cube.esti_error;
            split_direction = source_cube.split_direction;
        }

    };

    template <typename Real>
    auto make_cube_2D(const Real& x_min, const Real& x_max, const Real& y_min, const Real& y_max){
        Cube<Real, 2> maked_cube(std::array<Real, 2>{0.5*(x_min+x_max), 0.5*(y_min+y_max)},  std::array<Real, 2>{(x_max-x_min), (y_max-y_min)});
        return maked_cube;
    }

    template <typename Real>
    auto make_cube_3D(const Real& x_min, const Real& x_max, const Real& y_min, const Real& y_max, const Real& z_min, const Real& z_max){
        Cube<Real, 3> maked_cube(std::array<Real, 3>{0.5*(x_min+x_max), 0.5*(y_min+y_max), 0.5*(z_min+z_max)},  \
            std::array<Real, 3>{(x_max-x_min), (y_max-y_min), (z_max-z_min)});
        return maked_cube;
    }

    template<typename Real, int N>
    auto split_cube(const Cube<Real, N>& cube){
        const std::array<Real, 2> split_policy{ 0.25, -0.25 };
        std::vector<Cube<Real, N> > small_cubes;

        if (cube.evaluated == false){
            throw std::invalid_argument("Current cube has not been evaluted yet.");
        }

        std::array<Real, N> cur_size = cube.cube_size;
        int j = cube.split_direction;
        // std::cout << "Split: " << cube.split_direction << std::endl;
        cur_size[j] /= 2.0;
        for (int i=0; i<2; i++){
            std::array<Real, N> cur_center = cube.cube_center;
            cur_center[j] += cube.cube_size[j] * split_policy[i];
            small_cubes.push_back( Cube<Real, N>(cur_center, cur_size) );
        }
        return small_cubes;
    }

    template<typename Real>
    auto divide_cube3D(const Cube<Real, 3>& cube){
        const std::array<std::array<Real, 3>, 8> divide_policy{
            std::array<Real, 3>{0.25, 0.25, 0.25},\
            std::array<Real, 3>{0.25, 0.25, -0.25},\
            std::array<Real, 3>{0.25, -0.25, 0.25},\
            std::array<Real, 3>{0.25, -0.25, -0.25},\
            std::array<Real, 3>{-0.25, 0.25, 0.25},\
            std::array<Real, 3>{-0.25, 0.25, -0.25},\
            std::array<Real, 3>{-0.25, -0.25, 0.25},\
            std::array<Real, 3>{-0.25, -0.25, -0.25}
        };
        
        if (cube.evaluated == false){
            throw std::invalid_argument("Current cube has not been evaluted yet.");
        }
        std::vector<Cube<Real, 3> > small_cubes;
        std::array<Real, 3> cur_size = cube.cube_size;
        for (int i=0; i<3; i++){
            cur_size[i] /= 2.0;
        }
        for (int i=0; i<8; i++){
            std::array<Real, 3> cur_center = cube.cube_center;
            cur_center[0] += cube.cube_size[0] * divide_policy.at(i).at(0);
            cur_center[1] += cube.cube_size[1] * divide_policy.at(i).at(1);
            cur_center[2] += cube.cube_size[2] * divide_policy.at(i).at(2);
            small_cubes.push_back( Cube<Real, 3>(cur_center, cur_size) );
        }
        return small_cubes;
    }

};