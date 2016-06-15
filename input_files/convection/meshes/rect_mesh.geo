// Gmsh project created on Thu Jun 11 14:07:28 2015

// The 4 points created

Point(1) = {0, 0, 0};
Point(2) = {1000, 0, 0};
Point(3) = {1000, 500, 0};
Point(4) = {0, 500, 0};

nb_cells_x = 20;
nb_cells_y = 10;

// Lines connecting the 4 points

Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

Transfinite Line{1,3} = nb_cells_x + 1;
Transfinite Line{2,4} = nb_cells_y + 1;

  
// Line loop for lines 1-4

Line Loop(6) = {1, 2, 3, 4};

// Ruled surface made from line loop 6. *Note that ruled surface defaults a square mesh and the plane surface defaults a triangular mesh.

Ruled Surface(6) = {6};

Transfinite Surface{6} = {1,2,3,4}; // points indices, ordered
Recombine Surface {6};

// Adding physical lines

Physical Line("top") = {1};
Physical Line("left") = {2};
Physical Line("bottom") = {3};
Physical Line("right") = {4};
Physical Surface("rock") = {6};
