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

// RIS-UE link of the interference reflected by RISs of other cells (NLoS, pathloss exponent 4).
// This exponent reproduces simulation_alpha/output.json exactly; the description there lists it as a
// factor to vary (e.g. 3.5 for a partial line-of-sight link: distance^3 * sqrt(distance)).
inline double PL_dN(const double& x, const double& y, const double& theta){
    double distance = std::sqrt(y * y + x * x - 2.0 * x * y * std::cos(theta)) + 1.0;
    return 1.0 / (distance * distance * distance * distance);
}
