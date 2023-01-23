jdutc = 2449470.5;

jdut1 = jdutc2jdut1(jdutc, 0);

    Tut1 = ((jdut1 - 2451545.0)/36525);
    
    T = Tut1;

    lambda = 218.32  + 481267.8813*T + 6.29*sind(134.9+477198.85*T) - 1.27*sind(259.2-413335.38*T) + 0.66*sind(235.7 + 890534.23*T) + 0.21*sind(269.9 + 954397.7*T) - 0.19*sind(357.5 + 35999.05*T) - 0.11*sind(186.6 + 966404.05*T);

    phi = 5.13*sind(93.3 + 483202.03*T) + 0.28*sind(228.2 + 960400.87*T) - 0.28*sind(318.3 + 6003.18*T) - 0.17*sind(217.6 - 407332.20*T);

    p = 0.9508 + (0.0518*cosd(134.9 + 477198.85*T)) + (0.0095*cosd(259.2 - 413335.38*T)) + (0.0078*cosd(235.7 + 890534.23*T)) + (0.0028*cosd(269.9 + 954397.70*T));

    ob = 23.439291 - 0.0130042*T - (1.64E-7)*T^2 + (5.04E-7)*T^3;

    rmoon = 1/sind(p); %output in AU

    rmod = rmoon*[cosd(phi)*cosd(lambda); cosd(ob)*cosd(phi)*sind(lambda) - sind(ob)*sind(phi); sind(ob)*cosd(phi)*sind(lambda) + cosd(ob)*sind(phi)];
    
zeta = 2306.2181*T + 0.30188*(T^2) + 0.017998*(T^3); % in arcseconds
theta = 2004.3109*T - 0.42665*(T^2) - 0.041833*(T^3);
z = 2306.2181*T + 1.09468*(T^2) + 0.018203*(T^3);

zeta = zeta*pi/648000; % change to radians from arcseconds
theta = theta*pi/648000;
z = z*pi/648000;

Qmodgcrf = R3(zeta)*R2(-theta)*R3(z);

rgcrf = Qmodgcrf*rmod;
posMoon = rgcrf;