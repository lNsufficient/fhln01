function K = getK(K_all, A, edof, nele, ndof)
%GETK getK(K_all, A, edof, nele, ndof)
%Calculates K matrix, when given K cell (for Area=1), Area of each element,
%edof matrix and number of elements and number of nodes.

K = zeros(ndof);
for i = 1:nele
    edof_ele = edof(i, 2:5);
    K(edof_ele, edof_ele) = K(edof_ele, edof_ele) + K_all{i}*A(i);
end

end

