#pragma once
#include <complex>


// API class: 
// * For practical implementation, it's sufficient to evaluate the Laplace 
//      transform's value at a specific point in the complex s-domain.
// ! However, virtual interfaces and templates are incompatible in C++, so 
//      I've implemented both complex and real value cases and will use a wrapper to bridge them.

// TODO: To enhance the robustness of the software, we may consider implementing mechanisms 
//      to address the following theoretical constraints:
// 1, Analytically verifying that the input values reside within the region of convergence 
//      of the Laplace transform could ensure meaningful results.
// 2, For extremely large input values, the `std::exp(x)` function may become numerically 
//      unstable. Implementing alternative approaches or limiting the input range might mitigate this issue.



class LaplaceTransform {
    public:
    virtual double eval(const double& s) const = 0;
    virtual std::complex<double> eval(const std::complex<double>& s) const = 0;

    // For computing derivatives
    
    // For computing exponents
    virtual double eval_exponent(const double& s) const = 0;
    virtual std::complex<double> eval_exponent(const std::complex<double>& s) const = 0;
};