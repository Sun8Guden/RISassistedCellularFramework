#pragma once

namespace risgeo {



// Global, compile-time configuration of the simulation (system constants and numerical truncations).

class Parameter{

    public:
        // Some antenna or fading related parameters are fixed, so we don't have to calculate them repeatedly.
        static constexpr double transmission_power {1.0}; // 1 Watt
        static constexpr double antenna_gain {6.332573977646111e-05}; // For 3 GHz dipole antenna
        static constexpr double zeta_mean {0.821658900384983}; // For multiplication of two Rician with parameter 1
        static constexpr double zeta_variance {0.324876651418140};

        // In numerical systems, we do not use real infinity, but somehow very large value! One can change these values as their wish.
        static constexpr double infinity_cauchy_principal_value {40.0}; // Width of the integration range in log(z); exp(20) is large
        static constexpr double infinity_large_threshold {10000.0}; // We only consider the range of SINR thresholds up to 10000.0 with little accuracy loss
        static constexpr double infinity_large_distance {3000.0}; // We only consider the range of interference up to 3000.0 meters with little accuracy loss.
        static constexpr double guard_distance {50.0}; // Smallest serving-BS distance considered when averaging over the serving distance.

        static constexpr double expected_error {1e-6};
};

} // namespace risgeo
