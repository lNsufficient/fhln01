function [xstar, errors]=xstar(lambda,C, A_max, A_min)

    %xstar = C*lambda;
    xstar = sqrt(C/lambda);
    
    too_big = (xstar > A_max);
    xstar(too_big) = A_max;

    too_small = (xstar < A_min);
    xstar(too_small) = A_min;
    
    errors = too_big - too_small; %1 where too big, -1 where to small
    %there is no way it will be too big AND to small, so they won't cancel
    %out.

end
  


