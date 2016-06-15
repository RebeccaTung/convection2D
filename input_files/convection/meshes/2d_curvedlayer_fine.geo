// Gmsh project created on Tue Mar 22 09:44:33 2016

/***********************************************/
/* Parameters  */
scale = 1e2;
xmin = -2.5*scale;
xmax = 2.5*scale;
ymin = -1.5*scale;
y1 = -0.4*scale; // bottom of middle layer
y2 = 0.4*scale; // top of middle layer
ymax = 1.5*scale;

/****** MESH ***********/
nb_cells_x = 50; // nb cells in X dir along fault

nb_cells_y_top = 20;
nb_cells_y_middle = 12; // must be even
nb_cells_y_bottom = 20;

mesh_type = 2; // 2=quads (fully regular), 1=quads (non-regular), 0=triangles
// if mesh_type == 0
lc1 = 0.4*scale; // mesh characteristic length outside the fault (for triangular mesh)
lc2 = 0.02*scale; // mesh characteristic length in fault (for triangular mesh)
// if mesh_type == 1 or 2
progress_coeff = 0.92; // progression coefficient (for denser regular mesh towards fault)
// if mesh_type == 1
nb_cells_x2 = 4; // nb cells in X dir on the sides
/***********************************************/


Point(1) = {xmin, ymin, 0, lc1};
Point(2) = {xmax, ymin,  0, lc1};
Point(3) = {xmax, y1, 0, lc2};
Point(4) = {xmin, y1, 0, lc2};
Point(5) = {xmax, y2, 0, lc2};
Point(6) = {xmin, y2, 0, lc2};
Point(7) = {xmax, ymax, 0, lc1};
Point(8) = {xmin, ymax, 0, lc1};
Point(9) = {0, y2, 0, lc2};
Point(10) = {0, y1, 0, lc2};
Point(11) = {xmin, (y1+y2)/2., 0, lc2};
Point(12) = {xmax, (y1+y2)/2., 0, lc2};
Point(13) = {xmin/2, 0, 0, lc2};
Point(14) = {xmax/2, (y1)*2, 0, lc2};
Point(15) = {xmin/2, (y2)*2, 0, lc2};
Point(16) = {xmax/2, 0, 0, lc2};
Point(17) = {xmin/2, y2, 0, lc2};
Point(18) = {0, 0, 0, lc2};
Point(19) = {xmax/2, y1, 0, lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Spline(3) = {4,13,10,14,3};
Line(4) = {4,1};
Line(51) = {3,12};
Line(52) = {12,5};
Spline(6) = {6,15,9,16,5};
Line(71) = {6,11};
Line(72) = {11,4};
Line(8) = {5,7};
Line(9) = {7,8};
Line(10) = {8,6};
Spline(11) = {11,17,18,19,12};

Line Loop(1) = {1,2,-3,4};
Line Loop(21) = {71, 11, 52, -6};
Line Loop(22) = {3, 51, -11, 72};
Line Loop(3) = {9, 10, 6, 8};

If(mesh_type==0)
  Printf("Building non-regular mesh (with triangles)");
  Plane Surface(11) = {1};
  Plane Surface(121) = {21};
  Plane Surface(122) = {22};
  Plane Surface(13) = {3};
EndIf

If (mesh_type!=0)
  Printf("Buidling regular mesh (with quads)");
  Transfinite Line{2,-4} = nb_cells_y_bottom + 1 Using Progression progress_coeff;
  Transfinite Line{51,72} = nb_cells_y_middle/2 + 1;
  Transfinite Line{52,71} = nb_cells_y_middle/2 + 1;
  Transfinite Line{-8,10} = nb_cells_y_top + 1 Using Progression progress_coeff;

  If (mesh_type==2)
    Printf("Mesh fully regular");
    Transfinite Line{1,3,6,9,11} = nb_cells_x + 1;
  EndIf
  If (mesh_type==1)
    Printf("Mesh non-fully regular");
    Mesh.Smoothing = 4; // to get a 4 step Laplacian smoothing of the mesh
    Transfinite Line{3,6,11} = nb_cells_x + 1;
    Transfinite Line{1,9} = Floor(nb_cells_x2 + 1);
  EndIf

  Plane Surface(121) = {21};
  Transfinite Surface{121} = {12,5,6,11};
  Plane Surface(122) = {22};
  Transfinite Surface{122} = {3,12,11,4};

  If (mesh_type==1)
    Plane Surface(11) = {1};
    Plane Surface(13) = {3};
  EndIf
  If (mesh_type==2)
    Ruled Surface(11) = {1};
    Ruled Surface(13) = {3};
    Transfinite Surface{11} = {1,2,3,4}; // points indices, ordered
    Transfinite Surface{13} = {5,7,8,6};
  EndIf

  Recombine Surface {11};
  Recombine Surface {121};
  Recombine Surface {122};
  Recombine Surface {13};


EndIf

//Physical Line must start from 0
Physical Line(0) = {1}; // bottom
Physical Line(1) = {2,51,52,8}; // right hand side
Physical Line(2) = {9};  // top
Physical Line(3) = {10,71,72,4}; // left hand side
Physical Line(4) = {6}; // top contact line between layer and matrix
Physical Line(5) = {3}; // bottom contact line between layer and matrix

Physical Surface(0) = {11}; // bottom block
Physical Surface(1) = {121,122}; // middle block
Physical Surface(2) = {13}; // top block

Physical Point(6) = {11,12}; // middle points on LHS and RHS boundaries
Physical Point(7) = {7,8}; // top corners
Physical Point(8) = {1,2}; // bottom corners

Physical Line(10) = {71,72}; // LHS boundary middle layer
Physical Line(11) = {51,52}; // RHS boundary middle layer


Coherence;
Coherence;
Coherence;
