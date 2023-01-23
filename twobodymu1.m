function dY = twobodymu1(t, y)
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    dY(4:6) = (-1/(norm(y(1:3))^3))*y(1:3);
end
