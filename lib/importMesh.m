function BVP = importMesh(BVP)
% QUADSCHECK Processes and validates the quadrilateral mesh.
%
%   BVP = QUADSCHECK(BVP) reads, checks, and adjusts the mesh file for
%   a plate with a hole. It ensures valid Jacobian matrices for elements,
%   prepares surface data for visualization, and updates the BVP structure.
%
% Inputs:
%   BVP - Structure containing the mesh data (e.g., nodes and elements).
%
% Outputs:
%   BVP - Updated structure with validated mesh data.
%
% Notes:
%   - The mesh must consist of quadrilateral elements.
%   - Jacobian determinants are checked to ensure positive orientation.
%   - Surface coordinates are prepared for potential visualization.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Extract Mesh Data
% Load mesh data generated from Gmsh
Mesh_infPlate
% Node coordinates and element connectivity
nodeCoordinates = msh.POS(:, 1:end - 1); % Nodes (x, y, z excluded if present)
elementNodes = msh.QUADS(:, 1:end - 1); % Quadrilateral elements

% Orderings for Jacobian and element corrections
order_jacobi = [1, 4, 3, 2]; % Correction for negative Jacobian
order2D = [1, 2, 4, 3]; % Standard ordering for natural derivatives

% Natural derivatives for 2D quadrilateral elements
N1_d = [-0.25, 0.25, -0.25, 0.25]; % Derivative with respect to ξ
N2_d = [-0.25, -0.25, 0.25, 0.25]; % Derivative with respect to η
N1_d = N1_d(order2D);
N2_d = N2_d(order2D);
naturalDerivatives = [N1_d; N2_d]'; % Combine derivatives

%% Validate and Correct Mesh
numElements = size(elementNodes, 1);
correctedCount = 0; % Count of elements with corrected orientation

for e = 1:numElements
    indices = elementNodes(e, :); % Nodes for the current element
    
    % Compute the Jacobian for the current element
    [JacobianMatrix, ~, ~] = Jacobian2D(nodeCoordinates(indices, :), naturalDerivatives);
    detJ = det(JacobianMatrix); % Determinant of the Jacobian
    
    % Check if Jacobian determinant is negative
    if detJ < 0
      correctedCount = correctedCount + 1;
      elementNodes(e, :) = indices(order_jacobi); % Correct element orientation
    end
    
end

%% Prepare Data for Visualization
xSurf = makeSurf(elementNodes, nodeCoordinates(:, 1)); % Surface X-coordinates
ySurf = makeSurf(elementNodes, nodeCoordinates(:, 2)); % Surface Y-coordinates

%% Optional: Visualization (Commented Out)
% Uncomment for debugging or visualization
% figure;
% hh = patch(xSurf, ySurf, 'cyan');
% axis off; axis equal;
% hh.FaceColor = 'none';
% hh.EdgeColor = [0.3010, 0.7450, 0.9330];

%% Update BVP Structure
BVP.msh.nodeCoordinates = nodeCoordinates;
BVP.msh.elementNodes = elementNodes;
BVP.msh.xSurf = xSurf;
BVP.msh.ySurf = ySurf;

%% Display Summary
fprintf('Import mesh and validation completed. %d elements corrected.\n', correctedCount);
end
