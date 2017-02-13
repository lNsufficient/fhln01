clear;
addpath '../calfem-3.4/fem';

%% Setup problem parameters
E = 210*10^3*10^3; %Pa
A_max = 0.1; %m^2
maxMass = 1000; %kg
density = 7800; %kg/m^3
l = 10+sqrt(2)*4; %m total length of all bars

%% Set up geometry
geom

%% Plot the geometry
plotpar = [1 4 1];

%specify area for each bar. 

A0 = 1; %so that K =K_0*A
Area = ones(nelm, 1)*A0;
fac = 1; %area factor

myeldraw2(Ex, Ey, plotpar, Area, fac) %Draw the geometry

%% Set up the K matrix. 
K = zeros(ndof);

%Here, A = 1, so K = K0*A
Ep = [ones(nelm, 1)*E, ones(nelm, 1)];
K_all = cell(nelm,1);
for el = 1:nelm;
    Ke = bar2e(Ex(el, :), Ey(el, :), Ep(el,:));
    K_all{el} = Ke;
    K(Edof(el, 2:5), Edof(el, 2:5)) = K(Edof(el, 2:5), Edof(el, 2:5)) + Ke;
end
%% Set up F matrix
F = f_ext;

%% Solve the system for u
u = solveq(K,F,bc);
ed = extract(Edof, u);
magnfac = 10;
myeldisp2(Ex,Ey,ed,plotpar,magnfac,Area,fac)

%% Formulate optimization problem
A_min = 0;
A_max = 0.1;
V_max = maxMass/density;
A_init = maxMass/(l*density);

TOL = 1e-2;

%% Optimization:
while norm(x - x_old) > TOL
    u = solveq(K,F,bc);
    %% Calculate derivatives
    C = zeros(nelm,1);
    for i = 1:nelm
        edof = Edof(i, 2:5);
        le = Le(i);
        ue = u(edof);
        Ke0 = K_all{i};
        C(i) = le/(ue'*Ke0*ue);
    end
    
    
    
    %% Get the new x
    x_old = x;
    x = xstar(lambda);
    
    %% Calculate new K
    K = zeros(nelm);
    for i = 1:nelm
        edof = Edof(i, 2:5);
        K(edof, edof) = K(edof, edof) + K_all{i}*x(i);
    end
end


parameter1=1;
parameter2=2;
lambda_min=.1;
lambda_max=100;


lambdastar = fzero(@(lambda) dphidlambda(lambda,parameter1,parameter2),[lambda_min lambda_max])
