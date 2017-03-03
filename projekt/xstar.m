function [xstar, errors]=xstar(lambda,C, x, x_max, x_min, alpha)

    %xstar = C*lambda;
    eta = 1/(1+alpha);
    xstar = (C/lambda).^eta.*x;
    
    too_big = (xstar > x_max);
    xstar(too_big) = x_max;

    too_small = (xstar < x_min);
    xstar(too_small) = x_min;
    
    errors = too_big - too_small; %1 where too big, -1 where to small
    %there is no way it will be too big AND to small, so they won't cancel
    %out.

end
  


