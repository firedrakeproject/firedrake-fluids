dx = 200.0;
Point(1) = {0,0,0,dx};
Extrude {3e3,0,0} {
  Point{1}; Layers{3e3/dx};
}
Extrude {0,1e3,0} {
  Line{1}; Layers{1e3/dx};
}
// Reserve 1 and 2 for top and bottom of extruded mesh
// Outer ends of the channel (x=0 and x=3e3)
Physical Line(3) = {3};
Physical Line(4) = {4};
// Sides of the channel (y=0 and y=1e3)
Physical Line(5) = {1};
Physical Line(6) = {2};
// 
Physical Surface(1) = {5};
