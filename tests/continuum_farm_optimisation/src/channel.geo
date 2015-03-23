basin_x = 4e3;
basin_y = 4e3;
site_x = 1e3;
site_y = 1e3;
element_size = 200;
element_size_coarse = 200;

Point(1) = {0, 0, 0, element_size_coarse};
Point(2) = {basin_x, 0, 0, element_size_coarse};
Point(3) = {0, basin_y, 0, element_size_coarse};
Point(4) = {basin_x, basin_y, 0, element_size_coarse};

Point(5) = {(basin_x - site_x)/2, (basin_y - site_y)/2, 0, element_size};
Extrude{site_x, 0, 0} { Point{5}; Layers{site_x/element_size}; }
Extrude{0, site_y, 0} { Line{1}; Layers{site_y/element_size}; }

Line(6) = {1, 2};
Line(7) = {2, 4};
Line(8) = {4, 3};
Line(9) = {3, 1};
Line Loop(10) = {9, 6, 7, 8};
Line Loop(11) = {3, 2, -4, -1};
Plane Surface(12) = {10, 11};
Physical Surface(13) = {12, 5};

// Outflow
Physical Line(4) = {7};

// Inflow
Physical Line(3) = {9};

// Sides
Physical Line(5) = {8};
Physical Line(6) = {6};

