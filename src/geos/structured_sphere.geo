lc = 1.0;
r = 1.0;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {r, 0, 0, lc};
//+
Point(3) = {0, r, 0, lc};
//+
Point(4) = {-r, 0, 0, lc};
//+
Point(5) = {0, -r, 0, lc};
//+
Point(6) = {0, 0, r, lc};
//+
Point(7) = {0, 0, -r, lc};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {2, 1, 6};
//+
Circle(6) = {6, 1, 4};
//+
Circle(7) = {4, 1, 7};
//+
Circle(8) = {7, 1, 2};
//+
Circle(9) = {3, 1, 6};
//+
Circle(10) = {6, 1, 5};
//+
Circle(11) = {5, 1, 7};
//+
Circle(12) = {7, 1, 3};

//+
Line(13) = {6, 1};
//+
Line(14) = {5, 1};
//+
Line(15) = {4, 1};
//+
Line(16) = {3, 1};
//+
Line(17) = {2, 1};
//+
Line(18) = {7, 1};
//+
Curve Loop(1) = {14, -17, -4};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {17, -16, -1};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {16, -15, -2};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {15, -14, -3};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {15, -13, 6};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {15, -18, -7};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {18, -17, -8};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {17, -13, -5};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {13, -14, -10};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {14, -18, -11};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {18, -16, -12};
//+
Surface(11) = {11};
//+
Curve Loop(12) = {16, -13, -9};
//+
Surface(12) = {12};
//+
Curve Loop(13) = {10, 4, 5};
//+
Surface(13) = {13};
//+
Curve Loop(14) = {5, -9, -1};
//+
Surface(14) = {14};
//+
Curve Loop(15) = {9, 6, -2};
//+
Surface(15) = {15};
//+
Curve Loop(16) = {6, 3, -10};
//+
Surface(16) = {16};
//+
Surface(17) = {13};
//+
Curve Loop(17) = {11, 8, -4};
//+
Surface(18) = {17};
//+
Curve Loop(18) = {8, 1, -12};
//+
Surface(19) = {18};
//+
Curve Loop(19) = {12, 2, 7};
//+
Surface(20) = {19};
//+
Curve Loop(20) = {11, -7, 3};
//+
Surface(21) = {20};
//+
Surface Loop(1) = {8, 2, 12, 14};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {12, 5, 3, 15};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {5, 4, 9, 16};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {9, 8, 1, 13};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {7, 10, 1, 18};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {7, 2, 11, 19};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {11, 6, 3, 20};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {6, 10, 4, 21};
//+
Volume(8) = {8};
