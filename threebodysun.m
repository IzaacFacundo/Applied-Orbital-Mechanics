function dY = threebodysun(t, y)
    musun = 1.327E11;
    mue = 398600.4415;
    days = t / 86400;
    x = sun(days)*149597870.7; % AU to km
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    dY(4:6) = (-mue/(norm(y(1:3))^3))*y(1:3) + musun*(((x - y(1:3))/(norm(x - y(1:3)))^3) - (x/(norm(x))^3));
end
