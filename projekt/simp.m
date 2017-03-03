clear;
load MBBCoarseMesh

%% Plot the Geometry
figure(1);
clf;
plotpar = [1 4 1];
nen = 4;
[ex, ey] = coordxtr(edof,coord,dof,nen);
nele = size(edof,1);
ndof = length(F);


%specify area for each bar.
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry


%% Parameters

t = 0.05; %thickness chosen by us
V_box = 0.3*0.1*t;
V_max = 0.4*V_box;

x_max = 1;
x_min = 1e-3;

q = 1;
alpha = 2;

w = 0.3/42;%m, width, length over nbr elements in x direction
h = 0.1/14;%m, hidth, height over nbr elements in y direction
ae = t*w*h;


E = 210*10^3*10^6; %Pa %THis was chosen by us, but is okay. Steel.
nu = 0.27; %Steel, poissons ratio. 
%E = 210*10^3; %Pa %THis was chosen by us, but is okay. 
D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
%D = D(1:2,1:2);

%% Set up the K matrix. 
ep = [1, t];

x = 0.4*ones(nele, 1);

K_all = cell(nele,1);
for i = 1:nele
    K_all{i} = planqe(ex(i,:),ey(i,:),ep,D*1);

end

K = getK_sheet(K_all, x, q, edof, nele, ndof);
%

% K = zeros(ndof);
% 
% %Here, A = 1, so K = K0*A
% ep = [ones(nele, 1)*E, ones(nele, 1)];
% K_all = cell(nele,1);
% for el = 1:nele
%     Ke = flw2i4e(ex(e1,:),ey,ep,D)
%     K_all{el} = Ke;
%     %K(edof(el, 2:5), edof(el, 2:5)) = K(edof(el, 2:5), edof(el, 2:5)) + Ke;
% end
% 

%% Get the displacements

u = solveq(K, F, bc);
ed = extract(edof,u);
fac = 100000;

figure(2);
clf;

eldisp2(ex,ey,ed,plotpar,fac)

%% Set up the optimization problem:

TOL = 7e-6; %Try decreasing if there are problems.
%TOL = 1e-17; %Remember to compare this value to the current value of A_min.
%TOL = 1e-8; %Fewer elements hit A_min as loop breaks before.
tol_c = 1e-6; %This is only for debugging purposes.

lambda_min = 1e-9;
lambda_max = 1e9;

x_old = inf;
nbr_runs = 0;


%% Optimization
res = [];

while norm(x - x_old,2) > TOL
    
   %% Calculate the new K and corresponding u
   K = getK_sheet(K_all, x, q, edof, nele, ndof);
   
   u = solveq(K,F,bc);
   
   %% Calculate derivatives
   C = zeros(nele, 1);

   for i = 1:nele
       edof_ele = edof(i, 2:end);
       
       u_ele = u(edof_ele);
       Ke0 = K_all{i};
       C(i) = (u_ele'*Ke0*u_ele)/ae;
   end
   
   lambdastar = fzero(@(lambda) dphidlambda(lambda, ae*ones(nele,1), C, x, x_max, x_min, V_max, alpha),[lambda_min lambda_max])
    %% Get the new x
    x_old = x;
    [x, errors] = xstar(lambdastar, C, x, x_max, x_min, alpha);
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
end


%% Plot the displacements
x_opt = x;
K = getK_sheet(K_all, x, q, edof, nele, ndof);
u = solveq(K, F, bc);
ed = extract(edof,u);
fac = 100000;

figure(3);
clf;

eldisp2(ex,ey,ed,plotpar,fac)



