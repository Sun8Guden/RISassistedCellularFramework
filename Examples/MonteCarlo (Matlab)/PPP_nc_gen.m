function locations = PPP_nc_gen(PoissonNumber, Rmin, Rmax, center)
%PPP_NC_GEN This function generate PPP in an circular area.

    curNumPoints = poissrnd(PoissonNumber);
    locations = zeros(curNumPoints, 2);

    rho = sqrt(rand(curNumPoints,1)*(Rmax^2-Rmin^2) + Rmin^2);
    theta = 2*pi*(rand(curNumPoints,1));

    [locations(:,1), locations(:,2)] = pol2cart(theta, rho);

    locations(:,1) = locations(:,1) + center(1);
    locations(:,2) = locations(:,2) + center(2);
end

