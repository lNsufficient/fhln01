function dof = node_dof(i)
%NODE_DOF returnerar index i u-vektorn f�r motsvarande displacement (hos
%noden), allts� [x_index, y_index];

dof = [2*i-1 2*i];

end

