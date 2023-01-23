%% Problem 2

r = [100 5000 -6000];
mu = 398600.0000;
J2 = 0.001082626700000;
J3 = -2.532700000000000e-06;
radEarth = 6300.0000;

ap = j2j3accel(r, J2, J3, radEarth, mu)

%% Problem 3

rsat = [2781.000 5318.000 -5629.000];
rsun = [2379260.000 148079334.000 -1009936.000];
Amr = 0.1;
Cr = 1.5;

a = asrp(rsat, rsun, Cr, Amr)

%% Problem 4

h = 0:1:1000;
rho = zeros(1,1001);
h0 = zeros(1,1001);
H = zeros(1,1001);
for i = 1:1001
    [rho(i), h0(i), H(i)] = getDensityParams( h(i) );
end


density = rho.*exp(-(h-h0)./H);
semilogx(density,h)
xlabel('Density (kg/m^3)')
ylabel('Altitude (km)')
title('Altitude vs. Density')



%% Problem 5
clear
clc

a = 6800;
e = 0.005;
i = deg2rad(71);
raan = deg2rad(300);
argp = deg2rad(78);
ta = 0;
mu = 398600.0000;

[r,v] = oe2rv(mu,a,e,i,raan,ta,argp);

y0 = [r; v];
time = 2458200.5:(1/(3600*24)):(2458201.5);
odeoptions = odeset('RelTol', 1e-10,'AbsTol',1e-20);
[T,Y] = ode45(@fullprop,time,y0,odeoptions);

% Doesn't work. I need to debug the fullprop function further.



