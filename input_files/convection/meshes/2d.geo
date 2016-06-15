/***********************************************/
/* Parameters  */
xmin = 0;
xmax = 2.0;
ymin = -1;
ymax = 0;

/****** MESH ***********/
nb_cells_x = 40; // nb cells in X dir along fault
nb_cells_y = 30; //nb cells in Y direction

mesh_type = 2; // 2=quads (fully regular), 1=quads (non-regular), 0=triangles

/***********************************************/

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin,  0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Transfinite Line{1,3} = nb_cells_x + 1;
Transfinite Line{2,4} = nb_cells_y + 1;

Line Loop(11) = {1,2,3,4};

Ruled Surface (11) = {11};

Transfinite Surface{11} = {1,2,3,4}; // points indices, ordered
Recombine Surface {11};


//Physical Line must start from 0
Physical Line(0) = {1}; // bottom
Physical Line(1) = {2}; // right hand side
Physical Line(2) = {3};  // top
Physical Line(3) = {4}; // left hand side

Physical Surface(0) = {11}; //block

Physical Point(21) = {4}; //top left node corner
Physical Point(22) = {3}; //top right node corner