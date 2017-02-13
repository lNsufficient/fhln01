clear;
addpath '../calfem-3.4/fem';

%% Setup problem parameters
E = 210*10^3*10^3; %Pa
maxArea = 0.1; %m^2
maxMass = 1000; %kg
density = 7800; %kg/m^3
l = 10+sqrt(2)*4; %m total length of all bars

%% Set up geometry
geom

%% Plot the geometry
plotpar = [1 4 1];

%specify area for each bar. 
A = maxMass/(l*density);
Area = ones(nelm, 1)*A;
fac = 1; %area factor

myeldraw2(Ex, Ey, plotpar, Area, fac) %Draw the geometry

%% Set up the K matrix. 
K = zeros(ndof);

Ep = [ones(nelm, 1)*E, Area];

for el = 1:nelm;
    Ke = bar2e(Ex(el, :), Ey(el, :), Ep(el,:));
    K(Edof(el, 2:5), Edof(el, 2:5)) = K(Edof(el, 2:5), Edof(el, 2:5)) + Ke;
end

%% Solve the systm for u


parameter1=1;
parameter2=2;
lambda_min=.1;
lambda_max=100;


lambdastar = fzero(@(lambda) dphidlambda(lambda,parameter1,parameter2),[lambda_min lambda_max])
