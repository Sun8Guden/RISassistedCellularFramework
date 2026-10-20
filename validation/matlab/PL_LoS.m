function channelGain = PL_LoS(loc_Tx, loc_Rx)
%   PL_FUNC calculate the power of the channel gain 
%   Gain = AntennaGain/( (distance+1)^alpha )

    % For 3GHz, AntennaGain = ((3e8/3e9)/(4*pi))^2 = 6.3326e-05;
    channelGain = 6.3326e-05 * ((norm(loc_Tx-loc_Rx) + 1.0)...
        ^( - 3 )); % Here, LoS pathloss exponent is 3;
end

