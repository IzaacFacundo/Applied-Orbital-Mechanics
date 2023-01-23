function [ta, fa, ma] = propogateh(time_array, mu, hmag, e, tp)
    ma = zeros(1,length(time_array));
    fa = zeros(1,length(time_array));
    ta = zeros(1,length(time_array));
    for i = 1:length(time_array)
        ma(i) = (mu^2/hmag^3)*(e^2-1)^(3/2)*(time_array(i)-tp);
        guess = ma(i);
        ff = @ (ecc) (-ecc + e*sinh(ecc) - ma(i));
        syms ecc;
        dff = matlabFunction(diff(ff(ecc)));
        eccold = guess;
        
        for j = 1:20
                eccnew = eccold - ff(eccold)/dff(eccold);
                eccold = eccnew;
        end
        fa(i) = eccold;
        ta(i) = 2*atan2(sqrt(e+1)*tanh(fa(i)/2),sqrt(e-1));
        if ta(i) < 0
                ta(i) = 2*pi + ta(i);
        end
    end
end