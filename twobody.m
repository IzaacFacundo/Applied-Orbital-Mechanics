function dY = twobody(t, y)
    dY = zeros(6,1);
    dY(1:3) = y(4:6);
    dY(4:6) = (-398600.4415/(norm(y(1:3))^3))*y(1:3);
end

    
