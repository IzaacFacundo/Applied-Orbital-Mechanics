function [a,emag,i,raan,argp,ta] = rv2oe(r,v,mu) % This function assumes e < 1
    rmag = norm(r);
    vmag = norm(v);
    h = cross(r,v);
    hmag = norm(h);
    e = (1/mu)*((((vmag^2) - (mu/rmag))*r) -((dot(r,v))*v));
    emag = norm(e);
    if emag < 0.0001
        emag = 0;
    end
    p = (hmag^2)/mu;
    a = p/(1-emag^2);
    k = [0 0 1];
    hk = h(3);
    i = acos(hk/hmag);
    n = cross(k,h);
    nmag = norm(n);
    
    if emag <= 0 && i<= 0
        X = ('True anomaly, RAAN, and argument of periapsis are undefined!');
        disp(X);
        argp = 'undefined';
        raan = 'undefined';
        tl = acos(r(1)/rmag);
        if r(2) < 0
            tl = 2*pi - tl;
        end
        ta = tl;
        X = ['Element 6: (', num2str(ta), ' rads) is actually True Longitude.'];
        disp(X);
    elseif emag <=0
        X = ('True anomaly and argument of periapsis are undefined!');
        disp(X);
        raan = 2*pi - acos(n(1)/nmag);
        argp = 'undefined';
        arg = dot((n./nmag),r)/rmag;
        if abs(arg-1.0) < 0.00001 && arg > 1.0
            arg = 1;
        elseif abs(arg+1.0) < 0.00001 && arg < -1.0
            arg = -1;
        end
        u = acos(arg);
        ta = u;
        X = ['Element 6: (', num2str(ta), ' rads) is actually Argument of Latitude.'];
        disp(X);
    elseif i <= 0
        X = ('RAAN and argument of periapsis are undefined!');
        disp(X);
        raan = 'undefined';
        argp = 'undefined';
        lp = acos(e(1)/emag);
        if e(2) < 0
            lp = 2*pi - lp;
        end
        ta = lp;
        X = ['Element 6: (', num2str(ta), ' rads) is actually Longitude of Periapsis.'];
        disp(X);
    else
        raan = 2*pi - acos(n(1)/nmag);
        argp = acos(dot(n,e)/(nmag*emag));
        if e(3) < 0
            argp = 2*pi - argp;
        end
        ta = acos(dot(r,e)/(rmag*emag));
    end
end
        
