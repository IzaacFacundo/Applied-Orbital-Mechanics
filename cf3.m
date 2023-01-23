%% I will be using a Newton Raphson Solver for this homework


M1 = 2.6;% Mach Number
t = 15; % Theta in degrees
t = deg2rad(t);
g = 1.4; % gamma, the gas constant for air
R = 287; % ideal gas constant for air
cp = -(R*g)/(1-g); 

f = @(b)((2*cot(b)*((M1^2)*(sin(b))^2 - 1))/((M1^2)*(g + cos(2*b)) + 2) - tan(t)); % b is Beta
iter = 0;
syms bsym
df = matlabFunction(diff(f(bsym)));
bm_func = @(f,b)(b-f(b)/df(b));
err = 1;
tol = 0.001;
maxiter = 50;
start = 1;
bm = start;
bmold1 = bm;
bmold2 = bm+0.01;
while err>tol && iter<maxiter    
    iter = iter +1;
    bm = bm_func(f,bmold1);    
    err = abs((bm-bmold1)/bm)*100;    
    bmold2 = bmold1;    
    bmold1 = bm;
end

if t > 0
    b = rad2deg(bm); % b is beta in degrees
elseif t == 0
    b = 90;
    bm = pi/2;
end
M1n = M1*sin(bm); % bm is beta in radians

M2n = sqrt((1 + ((g-1)/2)*(M1n^2))/(g*(M1n^2) - ((g-1)/2)));
M2 = (M2n)/sin(bm-t);

p21 = 1 + (2*g*(M1n^2 - 1))/(g+1); % static pressure ratio
dels = cp*log(p21*(2 + (g-1)*M1n^2)/((g+1)*M1n^2)) - R*log(p21); % change in entropy
po21 = exp(-dels/R); % stagnation pressure ratio


