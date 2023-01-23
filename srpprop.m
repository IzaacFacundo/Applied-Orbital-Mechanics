function dY = srpprop(t,y)
    dY = zeros(6,1);
    days = t / 86400;
    sunPos = sun(days);
    satPos = y;
    re = 6378.1363;
    Cr = 1.5;
    Amr = 0.01;
    rs = norm(satPos(1:3));
    rsatsun = sunPos - satPos(1:3);
    if dot(satPos(1:3), sunPos/norm(sunPos)) < - sqrt(rs^2 - re^2)
        gamma = 0;
    else 
        gamma = 1;
    end
    dY(1:3) = y(4:6);
    dY(4:6) = (-398600.4415/(norm(satPos(1:3))^3))*satPos(1:3) + ((-1/1000) * gamma * (1367/299792458) * Cr * Amr * (rsatsun/norm(rsatsun)));
end