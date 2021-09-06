r = 1.0;
z = 1.0;
h = 0.2;
//+
Point(1) = {0, 0, -0.5*z, h};
//+
Point(2) = {r, 0, -0.5*z, h};
//+
Point(3) = {0, r, -0.5*z, h};
//+
Point(4) = {-r, 0, -0.5*z, h};
//+
Point(5) = {0, -r, -0.5*z, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Line(3) = {1, 4};
//+
Line(4) = {1, 5};
//+
Circle(5) = {2, 1, 3};
//+
Circle(6) = {3, 1, 4};
//+
Circle(7) = {4, 1, 5};
//+
Circle(8) = {5, 1, 2};
//+
Curve Loop(1) = {6, -3, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -5, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, -8, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 7, -4};
//+
Plane Surface(4) = {4};
//+
Extrude {0, 0, z} {
  Surface{1}; Surface{2}; Surface{4}; Surface{3}; Layers{3}; Recombine;
}
