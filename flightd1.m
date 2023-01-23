%% Flight D #2
Sw = 232;
cwbar = 7.04;
xwacle = 4.07;
xwle = 16.4;
xwac = xwacle + xwle;
iw = deg2rad(1);
CL0w = -0.0443;
CLaw = 5.08;
CMacw = -0.0175;
rho = 0.002048; %slug/cuft
v = 354.44; % feet per second
q = 0.5*rho*(v^2);



St = 54;
ctbar = 3.83;
xtacle = 2.79;
xtle = 36.9;
xtac = xtle + xtacle;
CL0t = 0;
CLat = 4.26;
CMact = 0;
Cldet = 1.8;

eta0 = deg2rad(0.642);
etaa = 0.426;
qtq = 0.9;

xcg = 10.56;

W = 9500;
% CMcg = 0;
% CMcg = CM0 + CMa*a + CMit*it + CMde*de;

xcgbar = xcg/cwbar;
xacwbar = xwac/cwbar;
xactbar = xtac/ctbar;

CM0 = CMacw + qtq*(St/Sw)*(ctbar/cwbar)*CMact + (xcgbar - xacwbar)*(CL0w+CLaw*iw)+qtq*(St/Sw)*(xcgbar-xactbar)*(CL0t-CLat*eta0);

CMa = (xcgbar - xacwbar)*CLaw + qtq*(St/Sw)*(xcgbar-xactbar)*CLat*(1-etaa);

CMit = qtq*(St/Sw)*(xcgbar-xactbar)*CLat;

CMde = qtq*(St/Sw)*(xcgbar-xactbar)*Cldet;

CLit = qtq*(St/Sw)*CLat;

CLa = CLaw + qtq*(St/Sw)*CLat*(1-etaa);

%CL0 = CL0w + CLaw*iw + qtq * (St/Sw) * (CL0t - CLat*eta0);
CL0 = 0.1;
A = [CLa CLit; CMa CMit];
gamma = deg2rad(-3);
B = [((W*cos(gamma))/(q*Sw))-CL0 (-CM0)]';
X = linsolve(A,B);
