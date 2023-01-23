function dY = atmdrag(t, y)
    radEarth = 6378.136300000;
    Cd = 2; % Change when the problem changes!
    Amr = 0.01; % Same! A/m ratio
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    omega = 2*pi/86164;
    h = norm(y(1:3)) - radEarth;
    [rho, h0, H] = getDensityParams( h );
    density = rho.*exp(-(h-h0)./H); %kg/m^3
    vrel = y(4:6) - cross([0 0 omega]',y(1:3));
    dY(4:6) = (-398600.4415/(norm(y(1:3))^3))*y(1:3) - 0.5*Cd*Amr*density*norm(vrel)*vrel*1000; 
end