function g1 = g1(l, x, V_max)
%G1 Evaluation of constraint g1

g1 = l'*x - V_max; %l alltid positiv, så det är bara att köra.

end

