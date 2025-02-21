function BVP = infinitePlate_FEM_postprocessLE(BVP)
% INFINITEPLATE_FEM_POSTPROCESSLE Post-processes the FEM solution for a linear elastic infinite plate.
%   This function computes stresses and strains at nodes and elements,
%   and interpolates them to the surface for visualization and analysis.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       proc.xiRange1, proc.xiRange2  - Gauss point ranges in xi1 and xi2 directions
%       msh.nElm                     - Number of elements
%       proc.LE.u                    - Displacement vector
%       msh.nodeCoordinates          - Node coordinates
%       msh.elementNodes, msh.elementDOFs - Element connectivity and DOFs
%       material.C_pstrain           - Material stiffness matrix (plane strain)
%       material.E                   - Elasticity matrix
%
% Outputs:
%   BVP - Updated BVP structure with computed stresses and strains
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Assign from BVP
xik1 = BVP.proc.xiRange1; % Gauss points in xi1 direction
xik2 = BVP.proc.xiRange2; % Gauss points in xi2 direction
numberElements = BVP.msh.nElm; % Number of elements
displacements = BVP.proc.LE.u; % Displacement vector
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
et = BVP.msh.elementNodes; % Element connectivity
ed = BVP.msh.elementDOFs; % Element DOFs
C = BVP.material.C_pstrain; % Material constitutive matrix
E = BVP.material.E; % Elastic modulus matrix

%% Initialize Stresses and Strains
SIG = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % Stress tensor
strain = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % Strain tensor

%% Loop over Elements
for e = 1:numberElements
    indice = et(e, :); % Element nodes
    elementDof = ed(e, :); % Element DOFs
    c = 1; % Counter for Gauss points
    
    for j = 1:size(xik2, 2)
        
        for i = 1:size(xik1, 2)
            % Compute shape function derivatives in natural coordinates
            N1_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 1);
            N2_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 2);
            naturalDerivatives = [N1_d; N2_d]';
            
            % Compute Jacobian and derivatives in physical coordinates
            [~, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
            
            % Construct B matrix for strain-displacement relationship
            B = zeros(3, size(elementDof, 2));
            B(1, 1:2:end) = XYDerivatives(:, 1)';
            B(2, 2:2:end) = XYDerivatives(:, 2)';
            B(3, 1:2:end) = XYDerivatives(:, 2)';
            B(3, 2:2:end) = XYDerivatives(:, 1)';
            
            % Compute strains and stresses at Gauss points
            EPS = B * displacements(elementDof, :); % Strain
            SIG(e, c, :) = E * C * EPS; % Stress
            strain(e, c, :) = EPS; % Store strain
            c = c + 1;
        end
        
    end
    
end

%% Compute Node-Averaged Values
uniqueNodes = unique(et); % Unique nodes for averaging
sigma_xx_tmp = SIG(:, :, 1); sigma_xx = zeros(size(uniqueNodes, 1), 1);
sigma_yy_tmp = SIG(:, :, 2); sigma_yy = zeros(size(uniqueNodes, 1), 1);
sigma_xy_tmp = SIG(:, :, 3); sigma_xy = zeros(size(uniqueNodes, 1), 1);
strain_xx_tmp = strain(:, :, 1); strain_xx = zeros(size(uniqueNodes, 1), 1);
strain_yy_tmp = strain(:, :, 2); strain_yy = zeros(size(uniqueNodes, 1), 1);
strain_xy_tmp = strain(:, :, 3); strain_xy = zeros(size(uniqueNodes, 1), 1);

% Loop through unique nodes and compute average values
for ii = 1:size(uniqueNodes, 1)
  sigma_xx(ii) = mean(sigma_xx_tmp(ismember(et, uniqueNodes(ii))));
  sigma_yy(ii) = mean(sigma_yy_tmp(ismember(et, uniqueNodes(ii))));
  sigma_xy(ii) = mean(sigma_xy_tmp(ismember(et, uniqueNodes(ii))));
  strain_xx(ii) = mean(strain_xx_tmp(ismember(et, uniqueNodes(ii))));
  strain_yy(ii) = mean(strain_yy_tmp(ismember(et, uniqueNodes(ii))));
  strain_xy(ii) = mean(strain_xy_tmp(ismember(et, uniqueNodes(ii))));
end

% Interpolate values to the surface for visualization
sigma_xx_Surf = makeSurf(et, sigma_xx);
sigma_yy_Surf = makeSurf(et, sigma_yy);
sigma_xy_Surf = makeSurf(et, sigma_xy);
strain_xx_Surf = makeSurf(et, strain_xx);
strain_yy_Surf = makeSurf(et, strain_yy);
strain_xy_Surf = makeSurf(et, strain_xy);

%% Assign Results Back to BVP
BVP.postproc.LE.sigma_xx = sigma_xx;
BVP.postproc.LE.sigma_yy = sigma_yy;
BVP.postproc.LE.sigma_xy = sigma_xy;
BVP.postproc.LE.sigma_xx_Surf = sigma_xx_Surf;
BVP.postproc.LE.sigma_yy_Surf = sigma_yy_Surf;
BVP.postproc.LE.sigma_xy_Surf = sigma_xy_Surf;
BVP.postproc.LE.strain_xx = strain_xx;
BVP.postproc.LE.strain_yy = strain_yy;
BVP.postproc.LE.strain_xy = strain_xy;
BVP.postproc.LE.strain_xx_Surf = strain_xx_Surf;
BVP.postproc.LE.strain_yy_Surf = strain_yy_Surf;
BVP.postproc.LE.strain_xy_Surf = strain_xy_Surf;

end
