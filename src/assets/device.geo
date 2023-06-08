// Gmsh project created on Tue Mar 29 17:18:18 2022OCC version   : 7.3.0

e = 1;
h = 5e-2;

Point(1) = {0, 0, 0, h};
Point(2) = {e, 0, 0, h};
Point(3) = {e, e, 0, h};
Point(4) = {0, e, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line("port1") = {1};
Physical Line("port2") = {3};
Physical Surface("device") = {1};