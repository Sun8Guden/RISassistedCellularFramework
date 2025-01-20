function signal = func_signal_BS(distance, UE_loc)

    % Direct link experience Relaigh fading and NLoS channel
    signal = exprnd(1) * PL_NLoS([distance, 0], UE_loc);

end

