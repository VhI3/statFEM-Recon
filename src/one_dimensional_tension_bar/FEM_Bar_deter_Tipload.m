function displacement = FEM_Bar_deter_Tipload(E, f_bar, A, nElm, elementNodes, nodeCoordinates, activeDOFs, GDof)
% FEM_BAR_DETER_TIPLOAD - Finite Element Method solution for a 1D tension bar under tip load.
%
% This function computes the displacement field for a one-dimensional tension bar
% subjected to a tip load using linear finite element analysis.
%
% Inputs:
%   E              - Young's modulus [MPa]
%   f_bar          - Magnitude of the applied tip load [N]
%   A              - Cross-sectional area [mm²]
%   nElm           - Number of finite elements
%   elementNodes   - (nElm x 2) Element connectivity matrix
%   nodeCoordinates - (nNodes x 1) Nodal coordinates [mm]
%   activeDOFs     - Active degrees of freedom indices
%   GDof           - Total degrees of freedom
%
% Outputs:
%   displacement   - (GDof x 1) Displacement vector [mm]
%
% Notes:
%   - The function uses linear Lagrange shape functions and two-point Gauss integration.
%   - Essential boundary conditions are applied at the fixed node.
%   - The stiffness matrix is assembled element-wise.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

% derivative of lagrangian shape functions
Nshape_derivate = @(xi) [-1/2 1/2];

fixedNode = 1; % Fixed node (Dirichlet BC)

stiffness = zeros(GDof, GDof); % Global stiffness matrix
force = zeros(GDof, 1); % Global force vector
[xi, w] = Gauss_int(2); % Gauss points and weights

% Assemble global stiffness matrix
for e = 1:nElm
    elementDof = elementNodes(e, :);
    ke = zeros(length(elementDof));
    for i = 1:length(xi)
        [JacobianMatrix, ~, B] = Jacobian1D(nodeCoordinates(elementDof), Nshape_derivate(xi(i)));
        ke = ke + w(i) * (B' * E * A * B) * det(JacobianMatrix);
    end
    stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
end

% Apply boundary conditions and loads
force(end) = f_bar; % Apply tip load at tip of the bar
Krr = stiffness(activeDOFs, activeDOFs);
Kru = stiffness(activeDOFs, fixedNode);
Rr = force(activeDOFs);

% Solve for displacements
displacement = zeros(GDof, 1);
Uu = displacement(fixedNode);
Ur = Krr \ (Rr - Kru * Uu);
displacement(activeDOFs) = Ur;
end
