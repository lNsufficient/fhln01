function K = getK_sheet(K_all, x, q, edof, nele, ndof)
%GETK_SHEET Gets the global K matrix for a plane sheet element, given the
%element matrices in K_all. q is the penalization factor. 
K = zeros(ndof, ndof); 

for i = 1:nele
    Ke = K_all{i}*x(i)^q;
    
    K(edof(i,2:end), edof(i,2:end)) = K(edof(i,2:end), edof(i,2:end)) + Ke;

end

end