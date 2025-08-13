Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {-1, 0, 0};
Point(5) = {0, -1, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Physical Curve("Circle") = {1,2,3,4};
Mesh.ElementOrder = 1;    // straight segments
// Mesh.SecondOrderLinear = 1; // (optional) if you ever go to 2nd order, keep edges straight
// Mesh.SaveAll = 0;       // default; only physical groups are saved

// Mesh 1; Save "circle_lines.msh2";