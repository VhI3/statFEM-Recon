// Define characteristic lengths
h1 = 0.1; // Fine mesh size for the circle
h2 = 2;   // Coarser mesh size for the larger geometry
r = 0.4;  // Radius of the circle
l = 4;    // Length of the square side

// Define points
Point(1) = {0, r, 0, h1};   // Top of the circle
Point(2) = {r, 0, 0, h1};   // Right of the circle
Point(3) = {l, 0, 0, h2};   // Bottom-right of the square
Point(4) = {l, l, 0, h2};   // Top-right of the square
Point(5) = {0, l, 0, h2};   // Top-left of the square
Point(6) = {0, 0, 0, h2};   // Bottom-left of the square (center of circle)

// Define geometry (circle and square)
Circle(1) = {1, 6, 2};  // Circle arc from Point 1 to Point 2, centered at Point 6
Line(2) = {2, 3};       // Line from Point 2 to Point 3
Line(3) = {3, 4};       // Line from Point 3 to Point 4
Line(4) = {4, 5};       // Line from Point 4 to Point 5
Line(5) = {5, 1};       // Line from Point 5 to Point 1

// Define physical groups (boundary conditions)
Physical Line("Boundary_u_1") = {3}; // Horizontal boundary (Line 3)
Physical Line("Boundary_u_0") = {5}; // Vertical boundary (Line 5)
Physical Line("Boundary_v_0") = {2}; // Other boundary (Line 2)

// Define surface and physical surface
Line Loop(1) = {1, 2, 3, 4, 5}; // Closed loop for the surface
Plane Surface(1) = {1};         // Define a plane surface from Line Loop 1
Physical Surface("Surface") = {1}; // Physical surface for mesh generation
