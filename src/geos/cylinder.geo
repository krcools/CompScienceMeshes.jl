r = 1.0;
z = 1.0;
h = 1.0;
//+
Point(1) = {0, 0, -z, h};
//+
Point(2) = {0, 0, z, h};
//+
Point(3) = {r, 0, -z, h};
//+
Point(4) = {r, 0, z, h};
//+
Point(5) = {0, r, -z, h};
//+
Point(6) = {0, r, z, h};
//+
Point(7) = {-r, 0, -z, h};
//+
Point(8) = {-r, 0, z, h};
//+
Point(9) = {0, -r, -z, h};
//+
Point(10) = {0, -r, z, h};
//+
Circle(1) = {9, 1, 3};
//+
Circle(2) = {3, 1, 5};
//+
Circle(3) = {5, 1, 7};
//+
Circle(4) = {7, 1, 9};
//+
Circle(5) = {4, 2, 6};
//+
Circle(6) = {6, 2, 8};
//+
Circle(7) = {8, 2, 10};
//+
Circle(8) = {10, 2, 4};
//+
Line(9) = {7, 8};
//+
Line(10) = {5, 6};
//+
Line(11) = {3, 4};
//+
Line(12) = {9, 10};
//+
Curve Loop(1) = {9, 7, -12, -4};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {12, 8, -11, -1};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {2, 10, -5, -11};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {3, 9, -6, -10};
//+
Surface(4) = {4};
//+
Line(13) = {8, 2};
//+
Line(14) = {2, 6};
//+
Line(15) = {2, 4};
//+
Line(16) = {2, 10};
//+
Line(17) = {1, 5};
//+
Line(18) = {1, 7};
//+
Line(19) = {1, 9};
//+
Line(20) = {1, 3};
//+
Line(21) = {1, 2};
//+
Curve Loop(5) = {13, 14, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {14, -5, -15};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {16, 8, -15};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {13, 16, -7};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {17, 3, -18};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {18, 4, -19};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {19, 1, -20};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {20, 2, -17};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {13, -21, 18, 9};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {15, -11, -20, 21};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {14, -10, -17, 21};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {21, 16, -12, -19};
//+
Plane Surface(16) = {16};
//+
Surface Loop(1) = {10, 16, 13, 8, 1};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {16, 11, 2, 7, 14};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {3, 12, 6, 14, 15};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {15, 9, 4, 5, 13};
//+
Volume(4) = {4};
