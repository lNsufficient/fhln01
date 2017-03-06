clear;

addpath('../calfem-3.4/fem/')

load_coarse = 1;
filter_case = 3;
start_case = 3;

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
nele = size(edof,1);
ndof = length(F);
nnod = max(max(enod(:,2:end))); %should be ndof/2;


%specify area for each bar.
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
title('Inital geometry')
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry


%% Parameters

t = 10*1e-3; %thickness chosen by us
V_box = 0.3*0.1*t;
V_max = 0.4*V_box;

x_max = 1;
x_min = 1e-4;

q = 3;
alpha = 2;

if load_coarse
    w = 0.3/42;%m, width, length over nbr elements in x direction
    h = 0.1/14;%m, hidth, height over nbr elements in y direction
    ae = t*w*h;
else
    h = 0.1/40;
    w = 0.3/120;
    ae = t*w*h;
end

E = 210*10^3*10^6; %Pa %THis was chosen by us, but is okay. Steel.
nu = 0.27; %Steel, poissons ratio. 
%E = 210*10^3; %Pa %THis was chosen by us, but is okay. 
D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
%D = D(1:2,1:2);

%% Set up the K matrix. 
ep = [1, t];



if filter_case == 3
    xdim = nnod;
else
    xdim = nele;
end
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
title('All elements equal density, deformed');
clf;

eldisp2(ex,ey,ed,plotpar,fac)

%% Set up the optimization problem:

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
    max_nbr_runs = 50;
end
%TOL = 1e-17; %Remember to compare this value to the current value of A_min.
%TOL = 1e-8; %Fewer elements hit A_min as loop breaks before.
tol_c = 1e-6; %This is only for debugging purposes.

lambda_min = 1e-9;
lambda_max = 1e5;

x_old = inf;
nbr_runs = 0;


%% Optimization
Res = [];





times = zeros(max_nbr_runs, 5);

%% Preparations for filtering
% Calculate mean coordinates for each element
coord_mean = [mean(ex,2), mean(ey,2)];
coord_mean = [coord_mean; mean(-ex,2), mean(ey,2)];


%% FILTER

if filter_case == 1
    M = eye(nele);
elseif filter_case == 2;
    M = zeros(nele, nele);
    R = w*2;
    for i = 1:nele
        coord_ele = coord_mean(i,:);
        %This could have been done more efficient, dist(a,b) = dist(b,a)...
        dists = sqrt((coord_mean(:,1) - coord_ele(1)).^2 + (coord_mean(:,2) - coord_ele(2)).^2);
        phis = 3/(pi*R^2)*max(0, 1-dists/R);
        phis = phis(1:nele) + phis(nele+1:end);
        M(i,:) = phis'/sum(phis);

        if mod(i,1e22) == 0
            figure(8);
            clf;
            fill(ex', ey', M(i,:))
            colorbar;
        end
    end 
    x = M*x;
elseif filter_case == 3
    r = 2*w;
    ep_filt = [1, 3];
    
    %Not sure if this has to be done elementwise or if the following is
    %okay
    
    K_filt = zeros(nnod, nnod);
    M_filt = zeros(nnod, nnod);
    T_filt = zeros(nnod, nele);
    q = 1;
    for i = 1:nele
        ele_ind = enod(i,2:end); %This should be correct - rho only affects these indices.
        [K_filt_ele, T_filt_ele] = flw2i4e(ex(i,:), ey(i,:), ep_filt, eye(2), q);
        T_filt(ele_ind,i) = T_filt_ele;
        K_filt(ele_ind, ele_ind) = K_filt(ele_ind, ele_ind) + K_filt_ele;
        M_filt_ele = flw2i4m(ex,ey,1); 
        M_filt(ele_ind, ele_ind) = M_filt(ele_ind, ele_ind) + M_filt_ele;
    end
    
    drhotildedrho = (K_filt+M_filt)\T_filt;
end

res = inf;


while res > TOL
    
%while nbr_runs < max_nbr_runs
   nbr_runs = nbr_runs + 1;
   %% Calculate the new K and corresponding u
   tic
   if filter_case == 3
        K = zeros(ndof, ndof);
        ep_K = [1,t,3];
        for i = 1:nele
            Ke = plani4e_rho(ex(i,:),ey(i,:),ep_K,D,x(enod(i,2:end))',q);
            edof_ele = edof(i, 2:end);
            K(edof_ele, edof_ele) = K(edof_ele, edof_ele) + Ke;
        end
   else
        K = getK_sheet(K_all, x, q, edof, nele, ndof);
   end
   times(nbr_runs,1) = toc;
   
   tic
   u = solveq(K,F,bc);
   times(nbr_runs,2) = toc;
   %% Calculate derivatives
   
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
           dgdr_ele = getdgdrhotilde_el(ex(i,:),ey(i,:),ep_dg,D,u(edof(i,2:end))',x(ele_ind)',q);
           dgdr(ele_ind) = dgdr(ele_ind) + dgdr_ele;
       end
       C = drhotildedrho'*dgdr;
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
%     if any(errors == 1)
%         disp('hits the upper limits....');
%     end
%     if any(errors == -1)
%         disp('hits the lower limits.....');
%     end

%     
%     disp(sprintf('Current run was: %d', nbr_runs));
    %% Filter the design variables
    x = M*x;


    res = norm(x-x_old,2);
    disp(res)
    Res = [Res; res];
end




%% Plot the displacements
%Filter the optimal x
x = M*x;

x_opt = x;
K = getK_sheet(K_all, x, q, edof, nele, ndof);
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

fill(ex', ey', x)
colorbar;


%% Plot the residual
figure(5);
clf;
subplot(2,1,1)
plot(Res);
subplot(2,1,2)
semilogy(Res)

%%
figure(6)
clf
hist(x);
save_str = sprintf('coarse_%d_TOL_%2.2g_alpha_%d', load_coarse, TOL, alpha)
savefig(sprintf('%s.fig',save_str))

%% Save optimal x
x_opt = x;

save('simp_x_opt', 'x_opt');
save(save_str, 'x_opt');
