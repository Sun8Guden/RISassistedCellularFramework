#pragma once
#include <complex>
#include <cmath>
#include <iostream>
#include <cassert>

namespace risgeo {

// To enhance performance in this computationally intensive application, a C++ guru
// recommended using a static method to initialize the ChiSquare object. This avoids 
// the overhead of repeated initialization, which can be significant when performed a 
// billion times. We welcome suggestions regarding the code.

// Please note that the sign of 's' depends on the type of signal: negative for
// interference and noise, positive for the signal.

// The *_D1 functions return the derivative dL/ds (as Exponential::eval_D1 does).

class ChiSquare {

public:
    template <typename ComplexType>
    static ComplexType eval_central(ComplexType s, double variance){
        ComplexType value = 1.0 / std::sqrt( 1.0 + 2.0 * variance * s );
        return value;
    }


    template <typename ComplexType>
    static ComplexType eval_central_D1(ComplexType s, double variance){
        ComplexType value = - variance / ( std::sqrt( 1.0 + 2.0 * variance * s ) * ( 1.0 + 2.0 * variance * s ) );
        return value;
    }

    template <typename ComplexType>
    static ComplexType eval_non_central(ComplexType s, double non_central_parameter, double variance){
        ComplexType value = std::exp( - s * non_central_parameter * variance /( 1.0 + 2.0 * variance * s) ) / std::sqrt( 1.0 + 2.0 * variance * s );
        return value;
    }

    template <typename ComplexType>
    static ComplexType eval_non_central_D1(ComplexType s, double non_central_parameter, double variance){
        ComplexType value = - variance * std::exp( - s * non_central_parameter * variance /( 1.0 + 2.0 * variance * s) ) * \
        (1.0 + non_central_parameter + 2.0 * variance * s) /\
        ( std::sqrt( 1.0 + 2.0 * variance * s ) * ( 1.0 + 2.0 * variance * s ) * ( 1.0 + 2.0 * variance * s ) );
        return value;
    }
};

} // namespace risgeo
