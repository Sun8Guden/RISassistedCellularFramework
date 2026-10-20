function signal = func_signal_RIS(BS_distance_to_US, UE_loc, RIS_cluster, num_elements_RIS)
    signal = 0;
    % RIS fading is the results of link level analysis.
    % In my work, I model the fading amplitude to be a non-central Gaussian, 
    % as the sum of the reflection of each element, and applying the 
    % central limit theorem. 

    % The mean (0.8217) and standard derivation (0.3249) of the random
    % fading of the reflected link, as the result of the product of two Rician 
    % faded link (with rician coefficient as 1) 
    chi_mu = 0.821658900384983 * num_elements_RIS;
    chi_std = sqrt(0.324876651418140 * num_elements_RIS);

    for ii = 1:size(RIS_cluster, 1)

        reflected_signal_power = PL_LoS([BS_distance_to_US, 0], RIS_cluster(ii, :)) * ...
            PL_LoS(RIS_cluster(ii, :), UE_loc); 
        
   % RIS fading (you can modify this random variable based on your research, feel free to play)
        RIS_fading = (chi_mu + randn(1) * chi_std)^2;

        signal = signal + RIS_fading * reflected_signal_power;

    end
end

