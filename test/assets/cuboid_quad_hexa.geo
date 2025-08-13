// -----------------------------------------------------------------------------
//  Rectangular box meshed with quadrangles (2D) and hexahedra (3D)
//  Adapted for transfinite + recombination
// -----------------------------------------------------------------------------

// Box dimensions
lx = 2.0; // length x direction
ly = 2.0; // length y direction
lz = 2.0; // length z direction

// Structured divisions (edit as needed)
nx = 10; // #cells along x
ny = 8;  // #cells along y
nz = 6;  // #cells along z

// ----- Points -----
// Upper points (z = +lz/2)
Point(1) = {+lx/2, +ly/2, +lz/2};
Point(2) = {-lx/2, +ly/2, +lz/2};
Point(3) = {-lx/2, -ly/2, +lz/2};
Point(4) = {+lx/2, -ly/2, +lz/2};

// Lower points (z = -lz/2)
Point(5) = {+lx/2, -ly/2, -lz/2};
Point(6) = {-lx/2, -ly/2, -lz/2};
Point(7) = {-lx/2, +ly/2, -lz/2};
Point(8) = {+lx/2, +ly/2, -lz/2};

// ----- Lines -----
// Top (z = +lz/2)
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};

// Bottom (z = -lz/2)
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 5};

// Vertical edges
Line(9)  = {4, 5};
Line(10) = {1, 8};
Line(11) = {2, 7};
Line(12) = {3, 6};

// ----- Curve loops / Surfaces -----
// Numbering as in a cube
Curve Loop(1) = {1, 2, 3, 4};            // Top
Curve Loop(2) = {-4, 9, -8, -10};        // Side
Curve Loop(3) = {12, -5, -9, -3};        // Side
Curve Loop(4) = {-7, -11, -1, 10};       // Side
Curve Loop(5) = {-6, -12, -2, 11};       // Side
Curve Loop(6) = {5, 6, 7, 8};            // Bottom

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// -----------------------------------------------------------------------------
// Transfinite (structured) + Recombine => quads on faces, hex in volume
// -----------------------------------------------------------------------------

// Set structured divisions per edge
// x-directed edges: 1,3,5,7
Transfinite Line {1, 3, 5, 7} = nx Using Progression 1;
// y-directed edges: 2,4,6,8
Transfinite Line {2, 4, 6, 8} = ny Using Progression 1;
// z-directed edges: 9,10,11,12
Transfinite Line {9, 10, 11, 12} = nz Using Progression 1;

// Make all faces/surfaces transfinite and recombine to quads
Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface  {1, 2, 3, 4, 5, 6};

// Structured hex volume
Transfinite Volume {1};
// (Optional; with transfinite this is usually enough, but harmless to keep)
Recombine Volume {1};

// -----------------------------------------------------------------------------
// Physical groups (unchanged)
// -----------------------------------------------------------------------------
Physical Surface("Top")    = {1};
Physical Surface("Bottom") = {6};
Physical Surface("Side")   = {2, 3, 4, 5};
Physical Volume("Interior") = {1};

// -----------------------------------------------------------------------------
// Mesh I/O hints (optional)
// -----------------------------------------------------------------------------
// Mesh.MshFileVersion = 2.2; // enable if you need legacy .msh2
// Mesh.RecombineAll = 1;     // global switch (not required since we specify above)

// You can also have Gmsh mesh on parse by uncommenting the next two lines:
// Mesh 3;
// Save "box_hex.msh";
