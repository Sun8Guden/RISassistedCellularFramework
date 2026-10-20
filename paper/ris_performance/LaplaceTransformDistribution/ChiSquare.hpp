#pragma once
#include <complex>
#include <cmath>
#include <iostream>
#include <cassert>
#include <type_traits>

// Pay attention to the sign of s:
// For interference and noise, this domain of convergence is defined as negetive.
// For signal, this domain is positive.  

class ChiSquare {

public:
    template <typename S_type>
    static S_type eval_central(S_type s, double variance){
        
        // if (std::is_same<S_type, std::complex<double>>::value){
        //     assert(s.real() * variance > -0.5 && "Laplace transform of ChiSquare is not defined within the region of convergence.");        
        // } 

        S_type value = 1.0 / std::sqrt( 1.0 + 2.0 * variance * s );
        return value;
        
    }


    template <typename S_type>
    static S_type eval_central_D1(S_type s, double variance){
        
        // if (std::is_same<S_type, std::complex<double>>::value){
        //     assert(s.real() * variance > -0.5 && "Laplace transform of ChiSquare is not defined within the region of convergence.");        
        // } 

        S_type value = variance / ( std::sqrt( 1.0 + 2.0 * variance * s ) * ( 1.0 + 2.0 * variance * s ) );
        return value;
        
    }

    template <typename S_type>
    static S_type eval_non_central(S_type s, double non_central_parameter, double variance){
 
        // if (std::is_same<S_type, std::complex<double>>::value){
        //     assert(s.real() * variance > -0.5 && "Laplace transform of ChiSquare is not defined within the region of convergence."); 
        // } 
        // auto val = - s * non_central_parameter /( 1.0 + 2.0 * variance * s);
        S_type value = std::exp( - s * non_central_parameter /( 1.0 + 2.0 * variance * s) )/std::sqrt( 1.0 + 2.0 * variance * s );
        return value;
    }

    template <typename S_type>
    static S_type eval_non_central_D1(S_type s, double non_central_parameter, double variance){
 
        // if (std::is_same<S_type, std::complex<double>>::value){
        //     assert(s.real() * variance > -0.5 && "Laplace transform of ChiSquare is not defined within the region of convergence."); 
        // } 
        // auto val = - s * non_central_parameter /( 1.0 + 2.0 * variance * s);
        S_type value = std::exp( - s * non_central_parameter /( 1.0 + 2.0 * variance * s) ) *\
        (non_central_parameter * (0.125 + 0.25 * variance * s) + variance * ( 0.125 + 0.5 * s * variance + 0.5 * s * s * variance * variance )) /\
        (std::sqrt( 1.0 + 2.0 * variance * s ) * ( 0.5 + variance * s ) * ( 0.5 + variance * s ) * ( 0.5 + variance * s ));
        return value;
    }
};