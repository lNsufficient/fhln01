function dphidlambda=dphidlambda(lambda, l, C, x, x_max, x_min, V_max, alpha)

    dphidlambda = g1(l, xstar(lambda, C, x, x_max, x_min, alpha), V_max);
end
