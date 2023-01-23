function posMoon = moon(jdutc)

    jdut1 = jdutc2jdut1(jdutc, 0);

    Tut1 = ((jdut1 - 2451545.0)/36525);
    
    T = Tut1;

    lambda = 218.32  + 481267.8813*T + 6.29*sind(134.9+477198.85*T) - 1.27*sind(259.2-413335.38*T) + 0.66*sind(235.7 + 890534.23*T) + 0.21*sind(269.9 + 954397.7*T) - 0.19*sind(357.5 + 35999.05*T) - 0.11*sind(186.6 + 966404.05*T);

    phi = 5.13*sind(93.3 + 483202.03*T) + 0.28*sind(228.2 + 960400.87*T) - 0.28*sind(318.3 + 6003.18*T) - 0.17*sind(217.6 - 407332.20*T);

    p = 0.9508 + (0.0518*cosd(134.9 + 477198.85*T)) + (0.0095*cosd(259.2 - 413335.38*T)) + (0.0078*cosd(235.7 + 890534.23*T)) + (0.0028*cosd(269.9 + 954397.70*T));

    ob = 23.439291 - 0.0130042*T - (1.64E-7)*T^2 + (5.04E-7)*T^3;

    rmoon = 1/sind(p); %output in earth radii

    rj2000 = rmoon*[cosd(phi)*cosd(lambda); cosd(ob)*cosd(phi)*sind(lambda) - sind(ob)*sind(phi); sind(ob)*cosd(phi)*sind(lambda) + cosd(ob)*sind(phi)];
    
    dela = 0.0146*pi/648000;
    e0 = -0.16617*pi/648000;
    n0 = -0.0068192*pi/648000;

    Q = R3(-dela)*R2(e0)*R1(n0);
    posMoon = Q*rj2000;
end
