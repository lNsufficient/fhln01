clear;

addpath('../calfem-3.4/fem/')

load_coarse = 1; %Defines which mesh to choose. 
filter_case = 2; %filter_case = 1 : no filter, filter_case = 2: density filter,...
% filter_case = 3: Helmhotlz filter
start_case = 1; %Defines which initial guess to choose

if load_coarse
    load MBBCoarseMesh
else
    load MBBFineMesh.mat
end

%% Plot the Geometry
figure(1);
clf;
plotpar = [1 4 1];
nen = 4;
[ex, ey] = coordxtr(edof,coord,dof,nen);
nele = size(edof,1); %Number of elements
ndof = length(F);    %Degrees of freedom
nnod = max(max(enod(:,2:end))); %should be ndof/2;


%specify area of element, just for plotting
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
title('Inital geometry')
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry


%% Parameters

t = 10*1e-3; %thickness chosen by us
V_box = 0.3*0.1*t; %Volume of full box 
V_max = 0.4*V_box;

x_max = 1; %Maximum density
x_min = 1e-4; %MInimal density

q = 3; %Penalization factor
alpha = 2; %Factor in the OC method

if load_coarse
    w = 0.3/42;%m, width, length over nbr elements in x direction
    h = 0.1/14;%m, hidth, height over nbr elements in y direction
    ae = t*w*h; %Area of an element
else
    h = 0.1/40; 
    w = 0.3/120;
    ae = t*w*h;
end

E = 210*10^3*10^6; %Pa %THis was chosen by us, Elastic modulus for steel.
nu = 0.27; %Steel, poissons ratio. 
D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]; %Constitutive matrix


%% Set up the K matrix. 
ep = [1, t];
xdim = nele;

% Different initial guesses: 
if start_case == 1
    x = 0.4*ones(xdim, 1);
elseif start_case == 2
    x = 0.1*ones(xdim, 1);
elseif start_case == 3
    x = 0.01*ones(xdim,1);
elseif start_case == 4
    load('simp_x_opt')
    x = x_opt; 
end
    
%Constructs element K-matrices.   
K_all = cell(nele,1);
for i = 1:nele
    K_all{i} = planqe(ex(i,:),ey(i,:),ep,D*1);

end

%Constructs global K matrix. 
K = getK_sheet(K_all, x, q, edof, nele, ndof);

%     Ke = flw2i4e(ex(e1,:),ey,ep,D) %Don't use this to set up the K
%     matrix. 


%% Get the displacements for unoptimized structure

u = solveq(K, F, bc);
ed = extract(edof,u);
fac = 100000;

figure(2);
title('All elements equal density, deformed');
clf;

eldisp2(ex,ey,ed,plotpar,fac)

%% Set up the optimization problem:
%Different tolerances needed in different cases:
if load_coarse
    if q == 1
        TOL = 1e-6; %Try decreasing if there are problems.
    elseif q == 3
        TOL = 1e-9;
    end
    if start_case == 4
        TOL = 1.6e-11;
    end
    max_nbr_runs = 300;
else
    TOL = 1.3e-1; %The tolerance is adjusted to number of elements
    TOL=1e-6;
    max_nbr_runs = 50; %For debugging
end
%TOL = 1e-17; %Remember to compare this value to the current value of A_min.
%TOL = 1e-8; %Fewer elements hit A_min as loop breaks before.

lambda_min = 1e-9;
lambda_max = 1e5;

x_old = inf;
nbr_runs = 0;

%% Optimization
Res = [];
G0 = F'*u; %Compliance
G1 = g1(ae*ones(nele, 1), x, V_max); %Constraint

times = zeros(max_nbr_runs, 5); %Just for code optimization

%% Preparations for filtering
% Calculate mean coordinates for each element
coord_mean = [mean(ex,2), mean(ey,2)];
coord_mean = [coord_mean; mean(-ex,2), mean(ey,2)];


%% FILTER

if filter_case == 1 %No filter
    M = eye(nele);
elseif filter_case == 2; %Density filter
    M = zeros(nele, nele);
    R = w*1.5; %Radius of the filter. 
    for i = 1:nele
        coord_ele = coord_mean(i,:);
        %This could have been done more efficient, dist(a,b) = dist(b,a)...
        dists = sqrt((coord_mean(:,1) - coord_ele(1)).^2 + (coord_mean(:,2) - coord_ele(2)).^2);
        phis = 3/(pi*R^2)*max(0, 1-dists/R); 
        phis = phis(1:nele) + phis(nele+1:end); %The mirror image also influence the filter
        M(i,:) = phis'/sum(phis); %Adjusts weigths
        
    end 
    M = sparse(M);
    x = M*x; %Filters the design variables
elseif filter_case == 3 %Helmholtz filter
    r = w*1.2; %Radius of the filter
    ep_filt = [1, 3];
    
    %Not sure if this has to be done elementwise or if the following is
    %okay
    
    K_filt = zeros(nnod, nnod);
    M_filt = zeros(nnod, nnod);
    T_filt = zeros(nnod, nele);

    for i = 1:nele
        ele_ind = enod(i,2:end); %This should be correct - rho only affects these indices.
        [K_filt_ele, T_filt_ele] = flw2i4e(ex(i,:), ey(i,:), ep_filt, eye(2), 1);
        T_filt(ele_ind,i) = T_filt_ele;
        K_filt(ele_ind, ele_ind) = K_filt(ele_ind, ele_ind) + K_filt_ele;
        M_filt_ele = flw2i4m(ex(i,:),ey(i,:),1); 
        M_filt(ele_ind, ele_ind) = M_filt(ele_ind, ele_ind) + M_filt_ele;
    end
    
    M = (K_filt*r^2+M_filt)\T_filt;
    rho_tilde = M*x;
end

res = inf;

%Optimization loop: 
while res > TOL
    
%while nbr_runs < max_nbr_runs
   nbr_runs = nbr_runs + 1;
   %% Calculate the new K and corresponding u
   tic
   %Calculates global K-matrix
   if filter_case == 3
        K = zeros(ndof, ndof);
        ep_K = [1,t,3];
        for i = 1:nele
            Ke = plani4e_rho(ex(i,:),ey(i,:),ep_K,D,rho_tilde(enod(i,2:end))',q);
            edof_ele = edof(i, 2:end);
            K(edof_ele, edof_ele) = K(edof_ele, edof_ele) + Ke;
        end
   else
        K = getK_sheet(K_all, x, q, edof, nele, ndof);
   end
   K = sparse(K);
   times(nbr_runs,1) = toc;
   
   tic
   u = solveq(K,F,bc);
   times(nbr_runs,2) = toc;
   %% Calculate sensitivities
   
   tic
   if (filter_case == 1 || filter_case ==2)
       C = zeros(nele, 1);
       for i = 1:nele
           edof_ele = edof(i, 2:end);

           u_ele = u(edof_ele);
           Ke0 = K_all{i};
           C(i) = (u_ele'*q*x(i)^(q-1)*Ke0*u_ele)/ae;
       end
       C = M'*C;
   elseif filter_case == 3
       dgdr = zeros(nnod, 1);
       ep_dg = [1, t,3];
       for i = 1:nele
           ele_ind = enod(i, 2:end);
           dgdr_ele = getdgdrhotilde_el(ex(i,:),ey(i,:),ep_dg,D,u(edof(i,2:end))',rho_tilde(ele_ind)',q);
           dgdr(ele_ind) = dgdr(ele_ind) + dgdr_ele;
       end
       C = M'*dgdr/ae;
   end
      
   times(nbr_runs,3) = toc;
   
   tic
   lambdastar = fzero(@(lambda) dphidlambda(lambda, ae*ones(nele,1), C, x, x_max, x_min, V_max, alpha),[lambda_min lambda_max])
   times(nbr_runs,4) = toc;
    %% Get the new x
    x_old = x;
    tic
    [x, errors] = xstar(lambdastar, C, x, x_max, x_min, alpha);
    times(nbr_runs,5) = toc;

    %% Filter the design variables
    if filter_case == 3
        rho_tilde = M*x;
    else
        x = M*x; %rho_tilde could have been used in all cases, we guess.
    end


    res = norm(x-x_old,2);
    disp(res)
    Res = [Res; res];
    G0 = [G0; F'*u];
    G1 = [G1; g1(ae*ones(nele, 1), x, V_max)];
end

if load_coarse
    str = 'coarse';
else
    str = 'fine';
end

savestr = sprintf('%s_filter_%d_q_%d_alpha_%d_res_%1.2g_nbrRuns_%d', str, filter_case, q, alpha, Res(end), nbr_runs);


%% Plot the displacements
%Filter the optimal x

%THE FOLLOWING LINE MIGHT BE NEEDED IN filter_case 1 AND 2.
%x = M*x;

x_opt = x;
if filter_case == 3
    K = K; %WE SHOULD RECALCULATE THE K MATRIX HERE WHEN EVERYTHING WORKS!
else
    K = getK_sheet(K_all, x, q, edof, nele, ndof);
end
u = solveq(K, F, bc);
ed = extract(edof,u);
fac = 100000;

figure(3);

clf;

eldisp2(ex,ey,ed,plotpar,fac)
title(sprintf('Deformation of optimized structure,\n magnification factor %d',fac), 'interpreter', 'latex');

%%
figure(4);
clf;

fill([ex' -ex'], [ey' ey'], [x; x], 'linestyle', 'none')
colorbar;


%% Plot the residual
figure(5);
clf;
% subplot(2,1,1)
% plot(Res);
% subplot(2,1,2)
semilogy(Res)
xlabel('Iteration number');

%%
figure(6)
clf
hist(x);
save_str = sprintf('coarse_%d_TOL_%2.2g_alpha_%d', load_coarse, TOL, alpha)
savefig(sprintf('%s.fig',save_str))

%%
figure(7)

subplot(1, 2, 1)
plot(G0)
subplot(1, 2, 2)
plot(G1)

%% Save optimal x
x_opt = x;

save('simp_x_opt', 'x_opt');
save(save_str, 'x_opt');
