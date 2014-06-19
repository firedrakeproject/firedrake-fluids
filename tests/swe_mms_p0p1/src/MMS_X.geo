dx = 0.025;
Point(1) = {0,0,0,dx};
Extrude {1,0,0} {
  Point{1}; Layers{1/dx};
}
Extrude {0,1,0} {
  Line{1}; Layers{1/dx};
}

// Reserve 1 and 2 for top and bottom of extruded mesh

// x=0 and x=pi
Physical Line(3) = {3};
Physical Line(4) = {4};

// y=0 and y=pi
Physical Line(5) = {1};
Physical Line(6) = {2};

Physical Surface(1) = {5};
