function dY = fullprop(t,y)
    mu = 398600.4415000;
    J2 =  0.0010826267;
    J3 = -2.532700000000000e-06;
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    radEarth = 6378.136300000;
    j2j3acc = j2j3accel(y(1:3), J2, J3, radEarth, mu); % j2j3
    Amr = 0.01;
    Cr = 1.5;
    Cd =2;
    
    musun = 1.327E11;
    mumoon = 3903;
    
    days = t / 86400;
    x = sun(days)*149597870.7;
    asun = musun.*(((x - y(1:3))/(norm(x - y(1:3)))^3) - (x/(norm(x))^3));
    
    z = moon(days)*6378.1363; % earth radii to km
    amoon = mumoon.*(((z - y(1:3))/(norm(z - y(1:3)))^3) - (z/(norm(z))^3));


    omega = 7.2921e-5;
    h = norm(y(1:3)) - radEarth;
    [rho, h0, H] = getDensityParams( h );
    density = rho.*exp(-(h-h0)./H);
    vrel = y(4:6) - cross([0 0 omega]',y(1:3));
    ad = -0.5*Cd*Amr*density*norm(vrel)*vrel*1000;
    
    srp = asrp(y(1:3), x, Cr, Amr);
    
    dY(4) = (-mu/(norm(y(1:3))^3))*y(1) + j2j3acc(1) + ad(1) + srp(1) + amoon(1) + asun(1);
    dY(5) = (-mu/(norm(y(1:3))^3))*y(2) + j2j3acc(2) + ad(2) + srp(2) + amoon(2) + asun(2);
    dY(6) = (-mu/(norm(y(1:3))^3))*y(3) + j2j3acc(3) + ad(3) + srp(3) + amoon(3) + asun(3);
    
    
    
end
