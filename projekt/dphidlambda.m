function dphidlambda=dphidlambda(lambda, l, C, x, x_max, x_min, V_max, alpha)
% Calculates dphi/dlambda given a lambda, l  - the length of each element if
% bars are used - or the area of each element of each element if four node
% elements are used, the design variables x, x_max, x_min, V_max and alpha
% -  the parameter of the OC method. 
    dphidlambda = g1(l, xstar(lambda, C, x, x_max, x_min, alpha), V_max);
end
