function g1 = g1(l, x, V_max)
%G1 Evaluation of constraint g1

g1 = l'*x - V_max; %l alltid positiv, s� det �r bara att k�ra.

end

