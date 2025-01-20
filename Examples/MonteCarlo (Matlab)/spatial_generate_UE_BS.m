function [UE_loc, BS_locs] = spatial_generate_UE_BS(distance_UE, density_BS, area_radius)

    UE_loc = [0, 0];
    NumberBS = density_BS * pi * (area_radius^2 - distance_UE^2);    
    BS_locs = PPP_nc_gen(NumberBS, distance_UE, area_radius, [0, 0]);

end

