// Gmsh project created on Tue Jun  3 15:39:41 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.02};
//+
Point(2) = {2, 0, 0, 0.02};
//+
Point(3) = {1.5, 2, 0, 0.02};
//+
Point(4) = {0.5, 2, 0, 0.02};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
