// -----------------------------------------------------------------------------
//
//  This file is based on the Gmsh GEO tutorial 1 included in the Gmsh package
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------

// Define variables describing the dimensions of the inner and the outer conductor

lx = 2.0; // length x direction
ly = 2.0; // length y direction
lz = 2.0; // length z direction

// The variables can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point (which is omitted here, but can be used for local refinements):

// Upper points
Point(1) = {+lx/2, +ly/2, lz/2};
Point(2) = {-lx/2, +ly/2, lz/2};
Point(3) = {-lx/2, -ly/2, lz/2};
Point(4) = {+lx/2, -ly/2, lz/2};

// The distribution of the mesh element sizes will then be obtained by
// interpolation of these mesh sizes throughout the geometry. Another method to
// specify mesh sizes is to use general mesh size Fields (see `t10.geo'). A
// particular case is the use of a background mesh (see `t7.geo').

// If no target mesh size of provided, a default uniform coarse size will be
// used for the model, based on the overall model size.

// We can then define some additional points. All points should have different
// tags:

// Lower points
Point(5) = {+lx/2, -ly/2, -lz/2};
Point(6) = {-lx/2, -ly/2, -lz/2};
Point(7) = {-lx/2, +ly/2, -lz/2};
Point(8) = {+lx/2, +ly/2, -lz/2};

// Curves are Gmsh's second type of elementary entities, and, amongst curves,
// straight lines are the simplest. A straight line is identified by a tag and
// is defined by a list of two point tags. In the commands below, for example,
// the line 1 starts at point 1 and ends at point 2.
//
// Note that curve tags are separate from point tags - hence we can reuse tag
// `1' for our first curve. And as a general rule, elementary entity tags in
// Gmsh have to be unique per geometrical dimension.

// Top
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Bottom
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Bottom
Line(9) = {4, 5};
Line(10) = {1, 8};
Line(11) = {2, 7};
Line(12) = {3, 6};


// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has
// first to be defined. A curve loop is also identified by a tag (unique amongst
// curve loops) and defined by an ordered list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve to form
// a loop):

// Numbering as in a cube
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-4, 9, -8, -10};
Curve Loop(3) = {12, -5, -9, -3};
Curve Loop(4) = {-7, -11, -1, 10};
Curve Loop(5) = {-6, -12, -2, 11};
Curve Loop(6) = {5, 6, 7, 8};


// We can then define the surface as a list of curve loops

// As a general rule, if a surface has N holes, it is defined by N+1 curve loops:
// the first loop defines the exterior boundary; the other loops define the
// boundaries of the holes.

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set `Mesh.SaveAll=1;', or
// specify `-save_all' on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
//
// Here we define a physical curve for the inner and the outer conduct and
// and a physical surface with name "Interior" (with an automatic tag) 
// containing the geometrical surface 1:

Physical Surface("Top") = {1};
Physical Surface("Bottom") = {6};
Physical Surface("Side") = {2,3,4,5};

Physical Volume("Interior") = {1};

// Now that the geometry is complete, you can
// - either open this file with Gmsh and select `2D' in the `Mesh' module to
//   create a mesh; then select `Save' to save it to disk in the default format
//   (or use `File->Export' to export in other formats);
// - or run `gmsh capacitors_rectangular.geo -2 -format msh2` to mesh in 
// batch mode on the command line.

// You could also uncomment the following lines in this script:
//
//   Mesh 2;
//   Save "t1.msh"; // (S. Adrian: Use suffix .msh2 to get the msh2 format required by the code)
//
// which would lead Gmsh to mesh and save the mesh every time the file is
// parsed. (To simply parse the file from the command line, you can use `gmsh
// t1.geo -')

// By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
// format (the `MSH' format). You can save meshes in other mesh formats by
// specifying a filename with a different extension in the GUI, on the command
// line or in scripts. For example
//
//   Save "t1.unv";
//
// will save the mesh in the UNV format. You can also save the mesh in older
// versions of the MSH format:
//
// - In the GUI: open `File->Export', enter your `filename.msh' and then pick
//   the version in the dropdown menu.
// - On the command line: use the `-format' option (e.g. `gmsh file.geo -format
//   msh2 -2').
// - In a `.geo' script: add `Mesh.MshFileVersion = x.y;' for any version
//   number `x.y'.
// - As an alternative method, you can also not specify the format explicitly,
//   and just choose a filename with the `.msh2' or `.msh4' extension.
