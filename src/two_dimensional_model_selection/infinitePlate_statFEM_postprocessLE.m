function BVP = infinitePlate_statFEM_postprocessLE(BVP)
% INFINITEPLATE_STATFEM_POSTPROCESSLE Computes stress distribution postprocessing for statFEM.
%
% This function calculates stresses at integration points and averages them
% to nodes, providing stress distributions for visualization and analysis.
%
% Inputs:
%   BVP - Boundary value problem structure with statFEM results.
%
% Outputs:
%   BVP - Updated BVP structure containing postprocessed stress results.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
xik1 = BVP.proc.xiRange1; % Integration points in xi1 direction
xik2 = BVP.proc.xiRange2; % Integration points in xi2 direction
numberElements = BVP.msh.nElm; % Number of elements
displacements = BVP.statFEM.LE.mean_u_y; % Posterior mean displacements
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
et = BVP.msh.elementNodes; % Element connectivity
ed = BVP.msh.elementDOFs; % Element degrees of freedom
C = BVP.material.C_pstrain; % Material stiffness matrix
E = BVP.UQ.mu_E; % Mean Young's modulus

%% Calculation: Stresses at nodes
SIG = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % Stress tensor (xx, yy, xy)

for e = 1:numberElements
  indice = et(e, :); % Nodes of the element
  elementDof = ed(e, :); % DOFs of the element
  c = 1; % Counter for integration points
  
  for j = 1:size(xik2, 2)
    
    for i = 1:size(xik1, 2)
      % Compute derivatives of shape functions in natural coordinates
      N1_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 1);
      N2_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 2);
      naturalDerivatives = [N1_d; N2_d]';
      
      % Compute Jacobian matrix and derivatives in global coordinates
      [~, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
      
      % Assemble the strain-displacement (B) matrix
      B = zeros(3, size(elementDof, 2));
      B(1, 1:2:end) = XYDerivatives(:, 1)';
      B(2, 2:2:end) = XYDerivatives(:, 2)';
      B(3, 1:2:end) = XYDerivatives(:, 2)';
      B(3, 2:2:end) = XYDerivatives(:, 1)';
      
      % Compute strains and stresses
      EPS = B * displacements(elementDof, :); % Strain tensor
      SIG(e, c, :) = E * C * EPS; % Stress tensor
      c = c + 1; % Increment counter
    end
    
  end
  
end

% Average stresses on nodes
uniqueNodes = unique(et);
sigma_xx_tmp = SIG(:, :, 1); sigma_xx = zeros(size(uniqueNodes, 1), 1);
sigma_yy_tmp = SIG(:, :, 2); sigma_yy = zeros(size(uniqueNodes, 1), 1);
sigma_xy_tmp = SIG(:, :, 3); sigma_xy = zeros(size(uniqueNodes, 1), 1);

for ii = 1:size(uniqueNodes, 1)
  % Average stresses for the current node
  sigma_xx(ii) = mean(sigma_xx_tmp(ismember(et, uniqueNodes(ii))));
  sigma_yy(ii) = mean(sigma_yy_tmp(ismember(et, uniqueNodes(ii))));
  sigma_xy(ii) = mean(sigma_xy_tmp(ismember(et, uniqueNodes(ii))));
end

% Generate surface data for visualization
sigma_xx_Surf = makeSurf(et, sigma_xx);
sigma_yy_Surf = makeSurf(et, sigma_yy);
sigma_xy_Surf = makeSurf(et, sigma_xy);

%% Assign back to BVP
BVP.statFEM.postproc.LE.sigma_xx = sigma_xx;
BVP.statFEM.postproc.LE.sigma_yy = sigma_yy;
BVP.statFEM.postproc.LE.sigma_xy = sigma_xy;

BVP.statFEM.postproc.LE.sigma_xx_Surf = sigma_xx_Surf;
BVP.statFEM.postproc.LE.sigma_yy_Surf = sigma_yy_Surf;
BVP.statFEM.postproc.LE.sigma_xy_Surf = sigma_xy_Surf;

end
