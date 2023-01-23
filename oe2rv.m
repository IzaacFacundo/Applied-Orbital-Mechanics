
function [r,v] = oe2rv(mu,a,e,i,raan,ta,argp)
    
    p = a*(1-e^2);
    
    
    
    
  
    
    
    if e <= 0 && i <= 0
        r = p*[cos(argp) sin(argp) 0]';
        v = sqrt(mu/p)*[-sin(argp) cos(argp) 0]';
       
    elseif e <= 0
        rpqw = p*[cos(argp) sin(argp) 0]';
        vpqw = sqrt(mu/p)*[-sin(argp) cos(argp) 0]';
        
        Q = R3(-raan)*R1(-i);
        r = Q*rpqw;
        v = Q*vpqw;
       
    elseif i <= 0
        rpqw = p*[cos(ta) sin(ta) 0]';
        vpqw = sqrt(mu/p)*[-sin(ta) cos(ta) 0]';
        
        Q = R1(i)*R3(-argp);
        r = Q*rpqw;
        v = Q*vpqw;
      
    else
        
        rpqw = (p/(1+e*cos(ta)))*[cos(ta) sin(ta) 0]';
        vpqw = sqrt(mu/p)*[-sin(ta) e+cos(ta) 0]';
        
        Q = R3(-raan)*R1(-i)*R3(-argp);
        
        r = Q*rpqw;
        v = Q*vpqw;
    end
end
    