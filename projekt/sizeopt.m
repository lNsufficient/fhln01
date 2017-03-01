clear;

addpath('..\calfem-3.4\')
addpath('../calfem-3.4/')
load geomSO

%% Calculate element lengths:
le = sqrt((ex(:,1) - ex(:,2)).^2 + (ey(:,1) - ey(:,2)).^2);
l_tot = sum(le);

%% Plot the Geometry
figure(1);
clf;
plotpar = [1 4 1];

%specify area for each bar.
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry

%% Parameters

V_max = 2000*1e-9; %m^3 - supposed to be 2000 mm^3
d_max = 20*1e-3; %m, supposed to be 20 mm.
d_min = 0.01*1e-7; %m TEST VARYING THIS ONE!

E = 210*10^3*10^6; %Pa WARNING - THIS WAS GUESSED! 

%Calculate corresponding A_max and A_min respectively...
A_max = (d_max/2)^2*pi;
A_min = (d_min/2)^2*pi;

A_init = V_max/l_tot;



%A_max = 200;
%% Set up the K matrix. 
K = zeros(ndof);

%Here, A = 1, so K = K0*A
ep = [ones(nele, 1)*E, ones(nele, 1)];
K_all = cell(nele,1);
for el = 1:nele;
    Ke = bar2e(ex(el, :), ey(el, :), ep(el,:));
    K_all{el} = Ke;
    K(edof(el, 2:5), edof(el, 2:5)) = K(edof(el, 2:5), edof(el, 2:5)) + Ke;
end

x = ones(nele, 1)*A_init;
K = getK(K_all, x, edof, nele, ndof);

%% Set up F matrix
F = f;

%% Solve the system for u - NOT OPTIMIZED, just get an idea of what it looks like
u = solveq(K,F,bc);
ed = extract(edof, u);
magnfac = 1000;
figure(2);
clf;
myeldisp2(ex,ey,ed,plotpar,magnfac,x,fac)

%% Set up the optimization problem:

TOL = 1e-11; %Try decreasing if there are problems.
tol_c = 1e-6; %This is only for debugging purposes.

lambda_min = 1e-9;
lambda_max = 1e9;

x_old = inf;
nbr_runs = 0;

%% Optimization loop:

while norm(x - x_old,inf) > TOL
   %% Calculate the new K and corresponding u
   K = getK(K_all, x, edof, nele, ndof);
   
   u = solveq(K,F,bc);
   
   %% Calculate derivatives
   C = zeros(nele, 1);
   
   for i = 1:nele
       edof_ele = edof(i, 2:5);
       l_ele = le(i);
       u_ele = u(edof_ele);
       Ke0 = K_all{i};
       C(i) = (u_ele'*Ke0*u_ele)*x(i)^2/l_ele;
   end
   
%    if any(abs(C) < tol_c)
%        disp('C is very small');
%    end
   
    lambdastar = fzero(@(lambda) dphidlambda(lambda, le, C, A_max, A_min, V_max),[lambda_min lambda_max]);
    
    %% Get the new x
    x_old = x;
    [x, errors] = xstar(lambdastar, C, A_max, A_min);
    if any(errors == 1)
        disp('hits the upper limits....');
    end
    if any(errors == -1)
        disp('hits the lower limits.....');
    end
    nbr_runs = nbr_runs + 1;
    
    disp(sprintf('Current run was: %d', nbr_runs));
end

%% Solve the system for u
%K has to be updated with the new x.
K = getK(K_all, x, edof, nele, ndof);

u = solveq(K,F,bc);
ed = extract(edof, u);
magnfac = 1000;
plotpar = [1 4 1];
fac = 1000000;

figure(3);
clf;
myeldisp2(ex, ey, ed, plotpar, magnfac, x, fac); 

%% Calculate stresses
ep = [ones(nele, 1)*E, x];

es = zeros(size(x));
for i = 1:nele
    es(i) = bar2s(ex(i,:), ey(i,:), ep(i,:), ed(i,:));
end

sigma = es./x;

%% Calculate zero bars:
zero_ind = find(errors == -1);
edof_zero = edof(zero_ind, :);
ed = extract(edof_zero, u);

magnfac = 1000;
plotpar = [1 2 1];
fac = 1000000;

ex_zero = ex(zero_ind,:);
ey_zero = ey(zero_ind,:);


figure(3);
hold on;
myeldisp2(ex_zero, ey_zero, ed, plotpar, magnfac, x(zero_ind), fac); 
