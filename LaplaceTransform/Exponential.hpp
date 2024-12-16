#pragma once
#include <complex>
#include <cmath>
#include <iostream>
#include <cassert>

// Pay attention to the sign of s:
// For interference and noise, this domain of convergence is defined as negetive.
// For signal, this domain is positive.  

class Exponential {

public:
    template <typename ComplexType>
    static ComplexType eval(ComplexType s, double variance){
        ComplexType value = 1.0 / ( 1.0 + variance * s );
        return value;
    }

    template <typename ComplexType>
    static ComplexType eval_D1(ComplexType s, double variance){
        ComplexType value = - variance / ( ( 1.0 + variance * s ) * ( 1.0 + variance * s ) );
        return value;
    }
};