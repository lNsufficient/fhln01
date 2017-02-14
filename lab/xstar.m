function [xstar, errors]=xstar(lambda,C, A_max, A_min)

    %xstar = C*lambda;
    xstar = sqrt(C/lambda);
    
    too_big = (xstar > A_max);
    xstar(too_big) = A_max;

    too_small = (xstar < A_min);
    xstar(too_small) = A_min;
    
    errors = sum(too_big) + sum(too_small);

end
  


