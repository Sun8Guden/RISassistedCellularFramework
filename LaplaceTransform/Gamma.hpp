#pragma once
#include <complex>
#include <cmath>
#include <iostream>
#include <cassert>

// Pay attention to the sign of s:
// For interference and noise, this domain of convergence is defined as negetive.
// For signal, this domain is positive.  

class Gamma {

public:
    template <typename ComplexType>
    static ComplexType eval(ComplexType s, double variance, double exponent){
        ComplexType value = 1.0 / std::pow( 1.0 + variance * s, exponent );
        return value;
    }
};