F = 500; 
ndof = 18;

Dof = [(1:2:17)', (2:2:18)'];
nen = 2;

Coord = [0 0;
        1 0;
        2 0;
        3 0;
        4 0;
        3 1;
        2 1;
        1 1;
        0 1];


Enod = [1 1 2;
        2 2 3;
        3 3 4;
        4 4 5;
        5 5 6;
        6 6 7;
        7 7 8;
        8 8 9;
        9 4 6;
        10 4 7;
        11 3 7;
        12 3 8;
        13 2 8;
        14 2 9];

Edof = [Enod(:,1), node_dof(Enod(:,2)) node_dof(Enod(:,3))];

bc_nodes = [1, 9];
bc = [];
for i = bc_nodes
    bc = [bc; node_dof(i)', zeros(2,1)];
end

f_ext = zeros(ndof, 1);
f_ext(10) = -F;

[Ex, Ey] = coordxtr(Edof,Coord,Dof,nen);