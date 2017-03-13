clear;

addpath('..\calfem-3.4\')
addpath('../calfem-3.4/fem/')

%%
load geomSO

%% Calculate element lengths:
le = sqrt((ex(:,1) - ex(:,2)).^2 + (ey(:,1) - ey(:,2)).^2);
l_tot = sum(le);


%% Plot the Geometry
figure(1);
clf;
plotpar = [1 4 1];

%specify area for each bar just to illustrate initial geometry
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry

%% Parameters
V_max = 2000*1e-9; %m^3 - supposed to be 2000 mm^3
alpha = 1; %%to make sure that eta = 1/2 inside xstar, corresponding to a
%CONLIN approximation. 

%Restrictions on the bar diameter
d_max = 20*1e-3; %m, supposed to be 20 mm.
d_min = 0.01*1e-7; %m TEST VARYING THIS ONE!

E = 210*10^3*10^6; %Pa %This was chosen by us, this is the elastic modulus...
...for steel.

%Calculate corresponding A_max and A_min respectively...
A_max = (d_max/2)^2*pi;
A_min = (d_min/2)^2*pi;

%Initial guess
A_init = V_max/l_tot; %We are supposed to try a few different initial guesses here.
%A_init = 1e-8;

%% Set up the K matrix. 
K = zeros(ndof);

%Here, A = 1, so K = K0*A
ep = [ones(nele, 1)*E, ones(nele, 1)];
K_all = cell(nele,1);
%Element stifness matrices are created:
for el = 1:nele
    Ke = bar2e(ex(el, :), ey(el, :), ep(el,:));
    K_all{el} = Ke; 
end

x = ones(nele, 1)*A_init;
K = getK(K_all, x, edof, nele, ndof); %The global stifnees matrix is created

%% Set up F matrix
F = f; %Given from geometry file.

%% Solve the system for u - NOT OPTIMIZED, just get an idea of what it looks like
u = solveq(K,F,bc);
ed = extract(edof, u);
magnfac = 1000;
figure(2);
clf;
myeldisp2(ex,ey,ed,plotpar,magnfac,x,fac) %Illustrated unoptimized structure

%% Set up the optimization problem:

TOL = 1e-11; %Try varying if there are problems.
%TOL = 1e-8; %Fewer elements hit A_min as loop breaks before.

lambda_min = 1e-9;
lambda_max = 1e9;

x_old = inf;
nbr_runs = 0;

%% Optimization loop:
load_old_opt = 0;
if load_old_opt
    load 'current_best.mat';
    disp('===========USING OLD VALUES!!!!!===========');
    %changes x to old best x
    %changes nbr_runs to old nbr_runs
end

res = []; %Stores residual over time in the optimization loop.
G0 = F'*u;
%while norm(x - x_old,2) > TOL
while nbr_runs < 380
   %% Calculate the new K and corresponding u
   K = getK(K_all, x, edof, nele, ndof); 
   u = solveq(K,F,bc);
   
   %% Calculate sensitivities
   C = zeros(nele, 1); 
   xk = x;
   for i = 1:nele
       edof_ele = edof(i, 2:5);
       l_ele = le(i);
       u_ele = u(edof_ele);
       Ke0 = K_all{i};
       C(i) = (u_ele'*Ke0*u_ele)/l_ele; 
   end
   
   lambdastar = fzero(@(lambda) dphidlambda(lambda, le, C, xk, A_max, A_min, V_max, alpha),[lambda_min lambda_max]);
    
    %% Get the new x
    x_old = x;
    [x, errors] = xstar(lambdastar, C, xk, A_max, A_min, alpha);
%     if any(errors == 1)
%         disp('hits the upper limits....');
%     end
%     if any(errors == -1)
%         disp('hits the lower limits.....');
%     end
     nbr_runs = nbr_runs + 1;
%     
%     disp(sprintf('Current run was: %d', nbr_runs));
    res = [res; norm(x-x_old,2)];
    G0 =[G0; F'*u];
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
myeldisp2(ex, ey, ed, plotpar, magnfac, x, fac); %Plot of optimized structure.

%% Calculate stresses
ep = [ones(nele, 1)*E, x]; 

es = zeros(size(x));
for i = 1:nele
    es(i) = bar2s(ex(i,:), ey(i,:), ep(i,:), ed(i,:)); 
end

sigma = es./x;

%% Calculate bars where neither the the lower limit are reached. 
%The purpose is to visualieze bars with close to zero crossectional area. 
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

%% Investigates values for sigma
other_ind = find(errors == 0); %We find the indecis of elements in the correct interval for A, i.e.
%indeces of the elements that does neither reach A_max nor A_min.
deviation = norm(abs(sigma(other_ind))-sqrt(lambdastar*E), inf) %This should be zero. 
save('current_best', 'x', 'nbr_runs');

% g0_con = G0;
% %% Save g0
% paperPpt= [345 550/2];%550/2 because two pictures per page
% scalefactor = 0.49;
% figurePpt = round([100, 100 paperPpt*scalefactor]); %2 per col, 4 per row.
% gcf = figure(10)
% 
% load('g0_simp.mat')
% semilogy(G0,'b')
% hold on;
% semilogy(g0_con,'r')
% l = legend('OC, $\alpha = 2$, $q = 3$', 'CONLIN');
% set(l,'interpreter','latex');
% ylabel('Komplians','interpreter','latex')
% xlabel('Iteration number','interpreter','latex');
% set(gcf,'PaperUnits','points','PaperSize',figurePpt(3:4),'PaperPosition',figurePpt)
% str = 'g0_comp';
% print('-dpng',str)