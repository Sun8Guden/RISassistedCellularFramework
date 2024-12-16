#pragma once



// I have designed this as a singleton pattern to facilitate 
// the configuration of the simulation.

class Parameter{

    public:
        // Some antenna or fading related parameters are fixed, so we don't have to calculate then repreatedly.
        static constexpr double transmission_power {1.0}; // 1 Watt 
        static constexpr double antenna_gain {6.332573977646111e-05}; // For 3 GHz bipole antenna
        static constexpr double zeta_mean {0.821658900384983}; // For multiplication of two Rician with parameter 1
        static constexpr double zeta_variance {0.324876651418140};

        // In numerical systen , we do not use real infinity, but somehow very large value! One can change these values as their wish.
        static constexpr double infinity_cauchy_principal_value {40.0}; // Exponential of this value is large
        static constexpr double infinity_large_threshold {10000.0}; // We only consider the range of SNR up to 5000.0 with little accuracy loss
        static constexpr double infinity_large_distance {3000.0}; // We only consider the range of interference up to 3000.0 meters with little accuracy loss.
        static constexpr double guard_distance {50.0}; // We only consider the range of interference up to 3000.0 meters with little accuracy loss.

        static constexpr double expected_error {1e-6};
        static constexpr int max_number_function_call {100000}; 
};