function channelGain = PL_NLoS(loc_Tx, loc_Rx)
%   PL_FUNC calculate the channel gain provided by this fucntional.
%   Gain = AntennaGain/( (distance+1)^alpha )

    % For 3GHz, AntennaGain = ((3e8/3e9)/(4*pi))^2 = 6.3326e-05;
    channelGain = 6.3326e-05 * ((norm(loc_Tx-loc_Rx) + 1.0)...
        ^( - 4 )); % Here, NLoS pathloss exponent is 4;
end

