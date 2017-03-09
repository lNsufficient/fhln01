function [xstar, errors]=xstar(lambda,C, x, x_max, x_min, alpha)
    %Computes xstar given input lambda, the calculated sensitivities C,
    %the design variables x and alpha, the exponent in the OC-method. In CONLIN
    %alpha = 1.
   
    eta = 1/(1+alpha);
    xstar = (C/lambda).^eta.*x;
    
    too_big = (xstar > x_max);
    xstar(too_big) = x_max;

    too_small = (xstar < x_min);
    xstar(too_small) = x_min;
    
    errors = too_big - too_small; %1 where too big, -1 where to small
    %there is no way it will be too big AND to small, so they won't cancel
    %out. This vector is used in the main program to sort out which
    %elements that hit upper or lower limit.

end
  


