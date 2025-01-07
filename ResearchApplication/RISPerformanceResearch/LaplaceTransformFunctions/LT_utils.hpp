# pragma once
# include <cmath>


inline double PL_d(const double& x) {
    double distance = x + 1.0;
    return 1.0 / (distance * distance * distance * distance);
}

inline double PL_r1(const double& y){
    double distance = y + 1.0;
    return 1.0 / (distance * distance * distance);
}

inline double PL_r2(const double& x, const double& y, const double& theta){
    double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
    return 1.0 / (distance * distance * distance);
}

// inline double PL_d(const double& x) {
//     double distance_loss = 1.0/ (x * x * x * x);
//     return distance_loss > 1.0 ? 1.0 : distance_loss;
// }

// inline double PL_r1(const double& y){
//     double distance_loss = 1.0 / (y * y * y);
//     return distance_loss > 1.0 ? 1.0 : distance_loss;
// }

// inline double PL_r2(const double& x, const double& y, const double& theta){
//     double distance = std::sqrt( y * y + x * x - 2.0 * x * y * std::cos(theta) );
//     double distance_loss = 1.0 / (distance * distance * distance);
//     return distance_loss > 1.0 ? 1.0 : distance_loss;
// }