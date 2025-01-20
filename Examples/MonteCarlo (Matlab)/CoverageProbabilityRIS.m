monte_carlo_number = 1e5;
density_BS = 1e-5;
distance_UE = 100; % r=100 meters
coverage_threshold = 1;
noise_power = 1e-13; % Watt
area_max = 3000; % In simulation, if you put infinite area, you will generate 
                 % infinite points, which is impossible for Monte Carlo
                 % simulation. Therefore, we must simulate a sufficiently 
                 % large, yet finite, area. This area should be chosen so 
                 % that the interference from points beyond its boundaries 
                 % has a negligible impact on the simulated results
                 % Here, we select 3000>>100. 


%% Clustered Random RIS Deployment
RIS_per_Cluster = 5;
RIS_elements_panel = 200; 
cluster_radius_min = 10.0;
cluster_radius_max = 25.0;


%%


coverage_count = 0;
parfor i=1:monte_carlo_number

    % generate a set of randomly located interference
    [UE_loc, BS_locs] = spatial_generate_UE_BS(distance_UE, density_BS, area_max);
    RIS_locs = PPP_nc_gen(RIS_per_Cluster, cluster_radius_min, cluster_radius_max, UE_loc);
        

    signal_BS = func_signal_BS(distance_UE, UE_loc);
    signal_RIS = func_signal_RIS(distance_UE, UE_loc, RIS_locs, RIS_elements_panel);
    interference = func_interference_BS_UE(BS_locs, UE_loc);
    SINR = (signal_BS + signal_RIS) / (interference + noise_power);
    coverage_count = coverage_count + (SINR > coverage_threshold);
        
end


coverage_prob = coverage_count / monte_carlo_number;
