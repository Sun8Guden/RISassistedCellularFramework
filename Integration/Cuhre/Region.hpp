#pragma once
#include <vector>
#include <array>

template <int N>
class Region{
    
    public:
        std::array<double, N> m_center;
        std::array<double, N> m_width;
        
        bool evaluated { false };
        double estimated_error {0.0};
        int divide_direction {0};

        Region(){}
        Region(const std::array<double, N>& center_arg, const std::array<double, N>& size_arg):\
        m_center(center_arg), m_width(size_arg){};
        ~Region(){};


        std::vector<Region> split();
        
};


template<int N>
std::vector<Region<N> > Region<N>::split(){
    if (evaluated == false){
        throw "The region is not evaluated yet, it makes no sense to divide and conqure.";
    }
    std::vector<Region > sub_regions;

    m_width.at(divide_direction) /= 2.0;
    m_center.at(divide_direction) -= 0.5 * m_width.at(divide_direction);
    sub_regions.push_back( Region(m_center, m_width) ); // direct initialization.

    m_center.at(divide_direction) += m_width.at(divide_direction);
    sub_regions.push_back( Region(m_center, m_width) );
    return sub_regions;
}

