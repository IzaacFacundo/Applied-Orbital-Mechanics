function acc = j2j3accel(satPos, J2, J3, radEarth, mu)
    
   
    accj2 = [0 0 0];
    accj3 = [0 0 0];
    
    i = satPos(1);
    j = satPos(2);
    k = satPos(3);
    
    r = norm(satPos);
    
    accj2(1) = ((-1.5.*J2.*mu.*(radEarth.^2).*i)./(r.^5)).*(1 - 5.*((k.^2) ./ (r.^2)));   
    accj2(2) = ((-1.5*J2*mu*radEarth.^2*j)./(r.^5))*(1 - 5*(k.^2 / r.^2)); 
    accj2(3) = ((-1.5*J2*mu*radEarth.^2*k)./(r.^5))*(3 - 5*(k.^2 / r.^2)); 
    
    accj3(1) = ((-5*J3*mu*(radEarth^3)*i)./(2*r^7))*(3*k - ((7*k.^3)/r^2));
    accj3(2) = ((-5*J3*mu*(radEarth^3)*j)./(2*r^7))*(3*k - ((7*k.^3)/r^2));
    accj3(3) = ((-5*J3*mu*(radEarth^3))/(2*(r^7))) * ((6*(k.^2)) - (7*(k.^4)/(r^2)) - (.6*(r^2)));

    acc = accj2 + accj3;


end