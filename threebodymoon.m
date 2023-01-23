function dY = threebodymoon(t, y)
    mumoon = 3903;
    mue = 398600.4415;
    days = t / 86400;
    x = moon(days)*6378.1363; % Earth Radii to km
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    dY(4:6) = (-mue/(norm(y(1:3))^3))*y(1:3) + mumoon*(((x - y(1:3))/(norm(x - y(1:3)))^3) - (x/(norm(x))^3));
end