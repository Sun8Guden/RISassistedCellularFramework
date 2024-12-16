#pragma once
#include <complex>
#include <array>

class LaplaceTransform {
    public:
    virtual std::complex<double> lt_eval(const std::complex<double>& s, const double& threshold, const double& distance)= 0;
    // virtual std::array<std::complex<double>, 2> lt_eval_array(const std::complex<double>& s, const double& threshold, const double& distance) {
    //     std::array<std::complex<double>, 2> return_value;
    //     return return_value;
    // };
};