Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Point(6) = {1, 1, -0, 1.0};
Circle(1) = {5, 1, 2};
Delete {
  Point{3};
}
Delete {
  Point{4};
}
Delete {
  Point{6};
}
Point(6) = {-1, 0, -0, 1.0};
Point(7) = {0, -1, -0, 1.0};
Circle(2) = {5, 1, 6};
Circle(3) = {6, 1, 7};
Circle(4) = {7, 1, 2};
Point(8) = {1, 0, 5, 1.0};
Point(9) = {0, 0, 5, 1.0};
Point(10) = {0, 1, 5, 1.0};
Point(11) = {-1, 0, 5, 1.0};
Point(12) = {0, -1, 5, 1.0};
Line(5) = {10, 5};
Line(6) = {11, 6};
Delete {
  Line{6, 5};
}
Circle(5) = {11, 9, 10};
Circle(6) = {10, 9, 8};
Circle(7) = {8, 9, 12};
Circle(8) = {12, 9, 11};
Line Loop(9) = {5, 6, 7, 8};
Line Loop(10) = {2, 3, 4, -1};
Ruled Surface(11) = {9, 10};
Delete {
  Surface{11};
}
Delete {
  Point{9, 1};
}
Delete {
  Point{9, 1};
}
Line(11) = {8, 2};
Line(12) = {12, 7};
Line(13) = {11, 6};
Line(14) = {10, 5};
Line Loop(15) = {11, -1, -14, 6};
Ruled Surface(16) = {15};
Line Loop(17) = {14, 2, -13, 5};
Ruled Surface(18) = {17};
Line Loop(19) = {13, 3, -12, 8};
Ruled Surface(20) = {19};
Line Loop(21) = {12, 4, -11, 7};
Ruled Surface(22) = {21};
Plane Surface(23) = {9};
Plane Surface(24) = {10};
Surface Loop(25) = {24, 18, 16, 22, 20, 23};
Volume(26) = {25};
Physical Surface(27) = {24};
Physical Surface(28) = {23};
Physical Surface(29) = {20, 18, 16, 22};
