clear;
addpath '../calfem-3.4/fem';

%% Setup problem parameters
E = 210*10^3*10^6; %Pa
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

figure(1)
clf;
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
A_min = 1e-4;
A_max = 1e-1
V_max = maxMass/density;
A_init = maxMass/(l*density);

TOL = 1e-9;
tol_c = 1e-6;

%% Optimization:
lambda_min=1e-9;
lambda_max=1e9;

x = A_init*ones(nelm,1);


x_old = inf;
nbr_runs = 0;
while norm(x - x_old,inf) > TOL
    %% Calculate new K
    K = zeros(ndof);
    for i = 1:nelm
        edof = Edof(i, 2:5);
        K(edof, edof) = K(edof, edof) + K_all{i}*x(i);
    end
    
    u = solveq(K,F,bc);
    %% Calculate derivatives
    C = zeros(nelm,1);
    for i = 1:nelm
        edof = Edof(i, 2:5);
        le = Le(i);
        ue = u(edof);
        Ke0 = K_all{i};
        C(i) = (ue'*Ke0*ue)*x(i)^2/le;
    end
    if any(abs(C) < tol_c)
        disp('C is very small')
    end
    lambdastar = fzero(@(lambda) dphidlambda(lambda, Le, C, A_max, A_min, V_max),[0 lambda_max]);
    
    %% Get the new x
    x_old = x;
    [x, errors] = xstar(lambdastar,C, A_max, A_min);
    if errors
        disp('hits the limits....');
    end
    nbr_runs = nbr_runs+1;
end

%% Solve the system for u
%K has to be updated with the new x. 
K = zeros(ndof);
    for i = 1:nelm
        edof = Edof(i, 2:5);
        K(edof, edof) = K(edof, edof) + K_all{i}*x(i);
    end

u = solveq(K,F,bc);
Ed = extract(Edof, u);
magnfac = 1;
plotpar = [1 3 1];
fac = 1000;

figure(2);
clf;
myeldisp2(Ex,Ey,Ed,plotpar,magnfac,x,fac)

Ep = [ones(nelm, 1)*E, x];

%DETTA FUNKAR INTE!!! :(
%es = bar2s(Ex,Ey,Ep,Ed);

%GÖR SÅ HÄR I STÄLLET!
Es = zeros(size(x));
for i = 1:nelm
    Es(i) = bar2s(Ex(i,:), Ey(i,:), Ep(i,:), Ed(i,:));
end

sigma = Es./x


