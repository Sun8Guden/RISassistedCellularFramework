# pragma once
# include <cmath>

// Since the pathloss function is known throughout the entire program, 
// it need not be explicitly passed between different objects.

// Future implementations may utilize alternative pathloss models to suit specific scenarios. 

class PathLoss {

    public:
    static double direct_nlos(const double& x) {
        double distance = x + 1.0;
        return 1.0 / (distance * distance * distance * distance);
    }

    static double partial_direct_nlos(const double& x) {
        double distance = x + 1.0;
        return - 4.0 / (distance * distance * distance * distance * distance);
    }

    static double direct_los(const double& x) {
        double distance = x + 1.0;
        return 1.0 / (distance * distance * distance);
    }

    static double reflected_los(const double& x, const double& y, const double& theta){
        double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
        return 1.0 / (distance * distance * distance); // 
    }

    static double reflected_los_derivative(const double& x, const double& y, const double& theta){
        double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
        double distance_minus_one = distance - 1.0;
        return - 3.0 * (x - y * std::cos(theta)) / (distance_minus_one * distance * distance * distance * distance); // 
    }


    // plos means partial line-of-sight, with pathloss exponent as 3.5
    static double reflected_plos(const double& x, const double& y, const double& theta){
        double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
        return 1.0 / (distance * distance * distance * std::sqrt(distance)); // 
    }

    static double reflected_nlos(const double& x, const double& y, const double& theta){
        double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
        return 1.0 / (distance * distance * distance * distance); // 
    }
};


// inline double PL_d(const double& x) {
//     double distance = x + 1.0;
//     return 1.0 / (distance * distance * distance * distance);
// }

// inline double PL_dN(const double& x, const double& y, const double& theta){
//     double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
//     return 1.0 / (distance * distance * distance * std::sqrt(distance)); // 
// }

// inline double PL_r1(const double& y){
//     double distance = y + 1.0;
//     return 1.0 / (distance * distance * distance);
// }

// inline double PL_r2(const double& x, const double& y, const double& theta){
//     double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
//     return 1.0 / (distance * distance * distance); // * std::sqrt(distance)
// }
