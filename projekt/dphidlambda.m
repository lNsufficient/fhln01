function dphidlambda=dphidlambda(lambda, l, C, A_max, A_min, V_max, alpha)

    dphidlambda = g1(l, xstar(lambda, C, A_max, A_min, alpha), V_max);
end
