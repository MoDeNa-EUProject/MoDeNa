// SetFactory("OpenCASCADE");

// Mesh.CharacteristicLengthMax = 0.4;

// R = 2;
// Block(1) = {0,0,0, R,R,R};
// s() = Abs(Boundary{Volume{1};});
// l() = Unique(Abs(Boundary{Surface{s()};}));
// pts() = Unique(Abs(Boundary{Line{l()};}));

// Characteristic Length{pts(0)} = 0.01;

// Periodic Surface{2} = {1} Translate{R,0,0};
// Periodic Surface{4} = {3} Translate{0,R,0};

// Physical Surface(1) = {1,2,3,4};
// Physical Surface(2) = {5,6};
// Physical Volume(1) = {1}

Point(1) = {-2,0,0};
Point(2) = {-1,0,0};
Point(3) = {-1,1,0};
Point(4) = {-2,1,0};
Point(5) = {-2,0,1};
Point(6) = {-1,0,1};
Point(7) = {-1,1,1};
Point(8) = {-2,1,1};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {4,8};
Line(11) = {2,6};
Line(12) = {3,7};
Line Loop(1) = {3,4,1,2};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {9,-8,-10,-4};
Line Loop(4) = {7,-10,1,12};
Line Loop(5) = {11,6,-12,2};
Line Loop(6) = {5,-11,3,9};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Surface Loop(1) = {1,2,3,-4,-5,-6};
Volume(1) = {1};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3,4,5,6};
Physical Volume(1) = {1};
