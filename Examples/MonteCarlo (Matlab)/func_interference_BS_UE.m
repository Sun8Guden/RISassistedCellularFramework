function Interference = func_interference_BS_UE(BS_locs, UE_loc)

    Interference = 0;
    for ii = 1:size(BS_locs, 1) % accumulate all the interference
        Interference = Interference + exprnd(1) * PL_NLoS(BS_locs(ii,:), UE_loc);
    end
end