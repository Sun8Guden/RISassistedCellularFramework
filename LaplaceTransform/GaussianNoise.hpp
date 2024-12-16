#pragma once
#include <cmath>
#include <complex>
// The Laplace transform of a constant is trivial, but I provide an object here.

class GaussianNoise {

public:
    template <typename ComplexType>
    static ComplexType eval(ComplexType s, double noise_power){
        ComplexType value = std::exp(- s * noise_power);
        // std::cout << "s is " << s << " noise is " << noise_power << " and value is " << value << std::endl;
        return value;
    }
};