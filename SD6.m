%% Problem 3.1
clear
clc

r = [-1.14697 0.75162 0.34193];
v = [0.65553 0.61048 0.44294];
mu = 1;

[a,emag,i,raan,argp,ta] = rv2oe(r,v,mu)

%% Problem 3.2
clear
clc

r = [21807.018 -18320.121 6800.183];
v = [-1.934079 -2.554491 -1.300623];
mu = 398600.4415;


[a,emag,i,raan,argp,ta] = rv2oe(r,v,mu)

%% Problem 3.3
clear
clc

mu = 398600.4415;
r = [21000*sqrt(2) 21000*sqrt(2) 0];
v = [-0.05*(sqrt(mu/210)) 0.05*sqrt(mu/210) 0];

[a,emag,i,raan,argp,ta] = rv2oe(r,v,mu)

%% Problem 4.1
clear
clc

mu = 398600.4415;
a = 34258.2;
e = 0.11112;
i = deg2rad(152.37);
raan = deg2rad(203.18);
argp = deg2rad(261.49);
ta = deg2rad(260.30);

[r,v] = oe2rv(mu,a,e,i,raan,ta,argp)

%% Problem 4.2
clear
clc

mu = 1;
a = 1.7;
e = 0.12645;
i = deg2rad(46.35);
raan = deg2rad(56.16);
argp = deg2rad(181.41);
ta = deg2rad(4.61);

[r,v] = oe2rv(mu,a,e,i,raan,ta,argp)



