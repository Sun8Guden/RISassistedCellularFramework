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
    template <typename S_type>
    static S_type eval(S_type s, double variance){
        // assert(s.real() * variance > -1.0 && "Laplace transform of Exponential is not defined within the region of convergence.");        

        S_type value = 1.0 / ( 1.0 + variance * s );
        return value;
    }

    template <typename S_type>
    static S_type eval_D1(S_type s, double variance){
        S_type value = variance /( ( 1.0 + variance * s ) * ( 1.0 + variance * s ) );
        return value;
    }

    template <typename S_type>
    static S_type eval_gamma(S_type s, double variance, double exponent){
        // assert(s.real() * variance > -1.0 && "Laplace transform of Gamma is not defined within the region of convergence.");        

        S_type value = 1.0 / std::pow( 1.0 + variance * s, exponent );
        return value;
    }
};