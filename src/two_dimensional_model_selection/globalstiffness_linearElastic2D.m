function [stiffness] = globalstiffness_linearElastic2D(BVP)
% GLOBALSTIFFNESS_LINEARELASTIC2D Computes the global stiffness matrix for a 2D linear elastic problem.
%   This function assembles the global stiffness matrix for a 2D finite
%   element method (FEM) problem using linear elasticity and plane strain.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       GDOFs            - Global degrees of freedom
%       msh.nElm         - Number of elements
%       msh.elementNodes - Element connectivity
%       msh.elementDOFs  - Element degrees of freedom
%       msh.nodeCoordinates - Nodal coordinates
%       proc.GP          - Gauss points and weights
%       material.C_pstrain - Material stiffness matrix for plane strain
%       geometry.T       - Thickness of the structure
%
% Outputs:
%   stiffness - Global stiffness matrix for the finite element problem
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Assign from BVP
GDOFs = BVP.msh.GDOFs; % Total number of global degrees of freedom
nElm = BVP.msh.nElm; % Number of elements
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
et = BVP.msh.elementNodes; % Element connectivity
ed = BVP.msh.elementDOFs; % Element degrees of freedom
xi1 = BVP.proc.GP.xi1; % Gauss points in xi1 direction
xi2 = BVP.proc.GP.xi2; % Gauss points in xi2 direction
w1 = BVP.proc.GP.w11; % Gauss weights in xi1 direction
w2 = BVP.proc.GP.w12; % Gauss weights in xi2 direction
thickness = BVP.geometry.T; % Plate thickness
C_pstrain = BVP.material.C_pstrain; % Material stiffness matrix (plane strain)

%% Initialization
stiffness = zeros(GDOFs, GDOFs); % Initialize global stiffness matrix

%% Loop over elements
for e = 1:nElm
  % Extract element connectivity and degrees of freedom
  indice = et(e, :); % Node indices for the element
  elementDof = ed(e, :); % Degrees of freedom for the element
  
  % Initialize element stiffness matrix
  ke = zeros(size(elementDof, 2));
  
  % Loop over Gauss points
  for j = 1:size(xi2, 2)
    
    for i = 1:size(xi1, 2)
      % Compute shape function derivatives in natural coordinates
      N1_d = Lagrange2D_d(xi1(i), xi2(j), 1, 1, [-1 1], [-1 1], 1);
      N2_d = Lagrange2D_d(xi1(i), xi2(j), 1, 1, [-1 1], [-1 1], 2);
      N1_d = N1_d([1 2 4 3]); % Reorder for consistent node numbering
      N2_d = N2_d([1 2 4 3]);
      naturalDerivatives = [N1_d; N2_d]';
      
      % Compute Jacobian and XY derivatives
      [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
      
      % Construct B matrix for strain-displacement relationship
      B = zeros(3, size(elementDof, 2));
      B(1, 1:2:end) = XYDerivatives(:, 1)';
      B(2, 2:2:end) = XYDerivatives(:, 2)';
      B(3, 1:2:end) = XYDerivatives(:, 2)';
      B(3, 2:2:end) = XYDerivatives(:, 1)';
      
      % Compute element stiffness matrix contribution
      ke = ke + w1(i) * w2(j) * (B' * C_pstrain * B) * det(JacobianMatrix) * thickness;
    end
    
  end
  
  % Assemble element stiffness matrix into global stiffness matrix
  stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
end

end
