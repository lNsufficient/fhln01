function dof = node_dof(i)
%NODE_DOF returnerar index i u-vektorn för motsvarande displacement (hos
%noden), alltså [x_index, y_index];

dof = [2*i-1 2*i];

end

