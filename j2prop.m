function dY = j2prop(t,satPos)
    mu = 398600.4415000;
    J2 =  0.0010826267;
    radEarth = 6378.136300000;
    dY = zeros(6,1);
    dY(1:3) = satPos(4:6);
    accj2 = [0 0 0];
    
    i = satPos(1);
    j = satPos(2);
    k = satPos(3);
    
    r = norm(satPos);
    
    accj2(1) = ((-1.5.*J2.*mu.*(radEarth.^2).*i)./(r.^5)).*(1 - 5.*((k.^2) ./ (r.^2)));   
    accj2(2) = ((-1.5*J2*mu*radEarth.^2*j)./(r.^5))*(1 - 5*(k.^2 / r.^2)); 
    accj2(3) = ((-1.5*J2*mu*radEarth.^2*k)./(r.^5))*(3 - 5*(k.^2 / r.^2)); 
   
    dY(4) = (-mu/(norm(satPos(1:3))^3))*i + accj2(1);
    dY(5) = (-mu/(norm(satPos(1:3))^3))*j + accj2(2);
    dY(6) = (-mu/(norm(satPos(1:3))^3))*k + accj2(3);

end