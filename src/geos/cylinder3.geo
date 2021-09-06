r = 1.0;
z = 1.0;
h = 1.0;
//+
Point(1) = {0, 0, -z, h};
//+
Point(2) = {r, 0, -z, h};
//+
Point(3) = {0, r, -z, h};
//+
Point(4) = {-r, 0, -z, h};
//+
Point(5) = {0, -r, -z, h};
//+
Point(6) = {r, 0, 0, h};
//+
Point(7) = {0, r, 0, h};
//+
Point(8) = {-r, 0, 0, h};
//+
Point(9) = {0, -r, 0, h};
//+
Point(10) = {0, 0, 0, h};
//+
Point(11) = {0, 0, z, h};
//+
Point(12) = {r, 0, z, h};
//+
Point(13) = {0, r, z, h};
//+
Point(14) = {-r, 0, z, h};
//+
Point(15) = {0, -r, z, h};
//+
Line(1) = {1, 5};
//+
Line(2) = {1, 4};
//+
Line(3) = {1, 3};
//+
Line(4) = {1, 2};
//+
Line(5) = {10, 7};
//+
Line(6) = {10, 9};
//+
Line(7) = {10, 6};
//+
Line(8) = {10, 8};
//+
Line(9) = {11, 13};
//+
Line(10) = {11, 12};
//+
Line(11) = {11, 15};
//+
Line(12) = {11, 14};
//+
Circle(13) = {4, 1, 5};
//+
Circle(14) = {5, 1, 2};
//+
Circle(15) = {2, 1, 3};
//+
Circle(16) = {3, 1, 4};
//+
Circle(17) = {9, 10, 6};
//+
Circle(18) = {6, 10, 7};
//+
Circle(19) = {7, 10, 8};
//+
Circle(20) = {8, 10, 9};
//+
Circle(21) = {13, 11, 12};
//+
Circle(22) = {12, 11, 15};
//+
Circle(23) = {15, 11, 14};
//+
Circle(24) = {14, 11, 13};
//+
Line(25) = {14, 8};
//+
Line(26) = {8, 4};
//+
Line(27) = {15, 9};
//+
Line(28) = {9, 5};
//+
Line(29) = {12, 6};
//+
Line(30) = {6, 2};
//+
Line(31) = {13, 7};
//+
Line(32) = {7, 3};
//+
Curve Loop(1) = {2, 13, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 14, -4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, 15, -3};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 16, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {6, 17, -7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, 18, -5};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {5, 19, -8};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {8, 20, -6};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {11, -22, -10};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {10, -21, -9};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {9, -24, -12};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {12, -23, -11};
//+
Plane Surface(12) = {12};
//+
Line(33) = {11, 10};
//+
Line(34) = {10, 1};
//+
Curve Loop(13) = {10, 29, -7, -33};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {33, 5, -31, -9};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {12, 25, -8, -33};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {27, -6, -33, 11};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {7, 30, -4, -34};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {34, 3, -32, -5};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {8, 26, -2, -34};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {6, 28, -1, -34};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {22, 27, 17, -29};
//+
Surface(21) = {21};
//+
Curve Loop(22) = {21, 29, 18, -31};
//+
Surface(22) = {22};
//+
Curve Loop(23) = {31, 19, -25, 24};
//+
Surface(23) = {23};
//+
Curve Loop(24) = {25, 20, -27, 23};
//+
Surface(24) = {24};
//+
Curve Loop(25) = {17, 30, -14, -28};
//+
Surface(25) = {25};
//+
Curve Loop(26) = {30, 15, -32, -18};
//+
Surface(26) = {26};
//+
Curve Loop(27) = {19, 26, -16, -32};
//+
Surface(27) = {27};
//+
Curve Loop(28) = {26, 13, -28, -20};
//+
Surface(28) = {28};
//+
Surface Loop(1) = {9, 21, 13, 16, 5};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {10, 22, 14, 6, 13};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {11, 23, 15, 14, 7};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {12, 24, 16, 15, 8};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {5, 2, 25, 17, 20};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {26, 3, 17, 6, 18};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {27, 4, 18, 7, 19};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {19, 28, 1, 20, 8};
//+
Volume(8) = {8};
