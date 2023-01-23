function  [ta, ea, ma] = propogate(time_array, n, e, tp)
    ma = zeros(1,length(time_array));
    ea = zeros(1,length(time_array));
    ta = zeros(1,length(time_array));
    if e == 0
       mes = ['The orbit is circular!'];
       disp(mes);
        for i = 1:size(time_array)
            ma(i) = n*(time_array(i)-tp);
            ea(i) = ma(i);
            ta(i) = ma(i);
        end
    elseif e < 1
        mes = ['The orbit is elliptical!'];
        disp(mes);
        
        for i = 1:length(time_array)
            ma(i) = n*(time_array(i)-tp);
            if ma(i) < pi
                guess = ma(i) + e/2;
            elseif ma(i) > pi
                guess = ma(i) - e/2;
            end
            fe = @ (ecc) (ecc - e*sin(ecc) - ma(i));
            syms ecc;
            dfe = matlabFunction(diff(fe(ecc)));
            eccold = guess;
            
            for j = 1:20
                eccnew = eccold - fe(eccold)/dfe(eccold);
                eccold = eccnew;
            end
            
            ea(i) = eccold;
            ta(i) = 2*atan2(sqrt(1+e)*tan(ea(i)/2),sqrt(1-e));
            if ta(i) < 0
                ta(i) = 2*pi + ta(i);
            end
        end
        
        
    elseif e > 1
        mes = ['The orbit is hyperbolic! Use propogateh.'];
        disp(mes);
    end
        end