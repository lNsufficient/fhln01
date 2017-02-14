function g0 = g0_conlin(F, x, u, K0, Edof)
%G0_CONLIN Evaluates g0 CONLIN approx

term1 = u'*K0*u;
for i = 1:size(Edof, 1)
    edof = Edof(i, 2:5);
    term1(edof, edof) = term1(edof, edof)*x(i); 
end




g0 = F'*u - 

end

