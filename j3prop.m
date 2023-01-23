function dY = j3prop(t,satPos)
    mu = 398600.4415000;
    J3 = -2.532700000000000e-06;
    dY = zeros(6,1);
    radEarth = 6378.136300000;
    dY(1:3) = satPos(4:6);
    accj3 = [0 0 0];
    
    i = satPos(1);
    j = satPos(2);
    k = satPos(3);
    
    r = norm(satPos);
    
    accj3(1) = ((-5*J3*mu*(radEarth^3)*i)./(2*r^7))*(3*k - ((7*k.^3)/r^2));
    accj3(2) = ((-5*J3*mu*(radEarth^3)*j)./(2*r^7))*(3*k - ((7*k.^3)/r^2));
    accj3(3) = ((-5*J3*mu*(radEarth^3))/(2*(r^7))) * ((6*(k.^2)) - (7*(k.^4)/(r^2)) - (.6*(r^2)));
   
    dY(4) = (-mu/(norm(satPos(1:3))^3))*i + accj3(1);
    dY(5) = (-mu/(norm(satPos(1:3))^3))*j + accj3(2);
    dY(6) = (-mu/(norm(satPos(1:3))^3))*k + accj3(3);

end