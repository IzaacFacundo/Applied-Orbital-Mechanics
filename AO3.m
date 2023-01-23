%% Problem 1.1, 1.2, 1.3
xp = 0.016768;
yp = 0.357565;
dut1 = -0.110914/86400;
dX = 0.000066;
dY = -0.000049;
X = 438.4313;
Y = 4.391260;
s = -0.006758;
timeutc = [02 19 2022 5 59 00];

% 1.1
jdut1 = gregutc2jdut1(timeutc, dut1); % 1.1
disp(jdut1);

% 1.2
tera = 2*pi*(0.7790572732640 + 1.00273781191135448*(jdut1 - 2451545.0)); 
disp(tera);

% 1.3
X = arcsec2deg(X);
Y = arcsec2deg(Y);
dX = arcsec2deg(dX);
dY = arcsec2deg(dY);
xp = arcsec2deg(xp);
yp = arcsec2deg(yp);
s = arcsec2deg(s);

    R = R3(-tera);
   disp(R);
   
   X = X+dX;
   Y = Y+dY;
   a = 0.5 + 0.125*(X^2 + Y^2);
   PN = [ 1-a*X*2 -a*X*Y X; -a*X*Y 1-a*Y^2 Y; -X -Y 1-a*(X^2+Y^2)];
   PN = deg2rad(PN)*R3(deg2rad(s));
   disp(PN);
   
   jdtt = (jdut1 - dut1) + (((69.184)/60)/60)/24;
   Ttt = ((jdtt - 2451545)/36525);
   sp = -0.000047*Ttt;
   W = R3(-sp)*R2(deg2rad(xp))*R1(deg2rad(yp));
   disp(W);
%% Problem 1.4
clear
clc

xp = 0.016768;
yp = 0.357565;
dut1 = -0.110914/86400;
dX = 0.000066;
dY = -0.000049;
X = 438.4313;
Y = 4.391260;
s = -0.006758;
timeutc = [02 19 2022 5 59 00];
jdut1 = gregutc2jdut1(timeutc,dut1);

Q = Qitrggcrf(jdut1, xp, yp, dX, dY, X, Y, s, dut1);
disp(Q); 


%% Problem 2
clear
clc

timeutc = [02 19 2022 5 59 00];
dut1 = -0.110914/86400;
jdut1 = gregutc2jdut1(timeutc,dut1);

% 2.1
ri = [-742.845 -5463.244 3196.066];
Tut1 = (jdut1 - 2451545) / 36525;
gmst = 67310.54841 + (876600*3600 + 8640184.812866)*Tut1 + 0.093104*(Tut1^2) - 0.0000062*(Tut1^3);
gmst = gmst/86400 - fix(gmst/86400);
gmstdeg = gmst*360;
disp(gmstdeg);

% 2.2
rg = R3(deg2rad(-gmstdeg))*ri';
disp(rg);

rggmst = rg;
xp = 0.016768;
yp = 0.357565;
dut1 = -0.110914/86400;
dX = 0.000066;
dY = -0.000049;
X = 438.4313;
Y = 4.391260;
s = -0.006758;
timeutc = [02 19 2022 5 59 00];
jdut1 = gregutc2jdut1(timeutc,dut1);

% 2.3
Q = Qitrggcrf(jdut1, xp, yp, dX, dY, X, Y, s, dut1);
rg = Q*ri';
disp(rg);

% 2.4
err = norm(rg) - norm(rggmst);
disp(err);
% The error between the two position vectors is massive. 
% It is over 400km away from where it should be. The GMST only
% transformation is unacceptable for our purposes.
