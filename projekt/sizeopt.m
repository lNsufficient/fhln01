clear;

addpath('..\calfem-3.4\')

load geomSO

%% Plot the Geometry
plotpar = [1 4 1];

%specify area for each bar.
A0 = 1;
Area = ones(nele, 1)*A0;
fac = 1; %area factor

figure(1)
clf;
myeldraw2(ex, ey, plotpar, Area, fac) %Draw the geometry

%%