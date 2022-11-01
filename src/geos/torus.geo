h=1.0;
//+
Point(1) = {4, 0, 0, h};
//+
Point(2) = {3, 0, 1, h};
//+
Point(3) = {3, 0, 0, h};
//+
Point(4) = {3, 0, -1, h};
//+
Point(5) = {2, 0, 0, h};
//+
Circle(1) = {2, 3, 1};
//+
Circle(2) = {1, 3, 4};
//+
Circle(3) = {4, 3, 5};
//+
Circle(4) = {5, 3, 2};

//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Curve{1}; Curve{2}; Curve{3}; Curve{4}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, 2*Pi/2} {
  Duplicata { Curve{1}; Curve{2}; Curve{3}; Curve{4}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, 3*Pi/2} {
  Duplicata { Curve{1}; Curve{2}; Curve{3}; Curve{4}; }
}
//+
Point(45) = {0, 0, 0, 1.0};
//+
Point(46) = {0, 0, 1, 1.0};
//+
Point(47) = {0, 0, -1, 1.0};
//+
Circle(17) = {5, 45, 18};
//+
Circle(18) = {18, 45, 31};
//+
Circle(19) = {31, 45, 44};
//+
Circle(20) = {44, 45, 5};
//+
Circle(21) = {13, 47, 26};
//+
Circle(22) = {26, 47, 39};
//+
Circle(23) = {39, 47, 4};
//+
Circle(24) = {4, 47, 13};
//+
Circle(25) = {6, 46, 19};
//+
Circle(26) = {19, 46, 32};
//+
Circle(27) = {32, 46, 2};
//+
Circle(28) = {2, 46, 6};
//+
Circle(29) = {21, 45, 34};
//+
Circle(30) = {34, 45, 1};
//+
Circle(31) = {1, 45, 8};
//+
Circle(32) = {8, 45, 21};
//+
Curve Loop(1) = {1, 31, -5, -28};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {5, 32, -9, -25};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {9, 29, -13, -26};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {13, 30, -1, -27};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {-17, -8, 28, 4};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {-18, -12, 25, 8};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {-19, -16, 26, 12};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {27, -4, -20, 16};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {-22, -15, 19, 11};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {-11, 18, 7, -21};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {-7, 17, 3, -24};
//+
Surface(11) = {11};
//+
Curve Loop(12) = {-3, 20, 15, -23};
//+
Surface(12) = {12};
//+
Curve Loop(13) = {-6, 24, 2, -31};
//+
Surface(13) = {13};
//+
Curve Loop(14) = {-30, -2, 23, 14};
//+
Surface(14) = {14};
//+
Curve Loop(15) = {-14, 22, 10, -29};
//+
Surface(15) = {15};
//+
Curve Loop(16) = {-10, 21, 6, -32};
//+
Surface(16) = {16};
