function BVP = infinitePlate_FEM_postprocessST(BVP)
% INFINITEPLATE_FEM_POSTPROCESSST Post-processes FEM results for a hyperelastic plate.
%   This function computes stresses and strains over time steps for a finite
%   element model of an infinite plate, including Cauchy stresses, PK1 stresses,
%   and strain components.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       Mesh, boundary conditions, material properties, and FEM results.
%
% Outputs:
%   BVP - Updated BVP structure with computed stresses and strains over time steps.
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Assign from BVP
xik1 = [-1 1]; % Gauss points in xi1 direction
xik2 = [-1 1]; % Gauss points in xi2 direction
numberElements = BVP.msh.nElm; % Number of elements
u_steps = BVP.proc.ST.u_steps; % Displacement steps
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
et = BVP.msh.elementNodes; % Element connectivity
ed = BVP.msh.elementDOFs; % Element degrees of freedom
E = BVP.material.E; % Young's modulus
C = BVP.material.C_pstrain; % Plane strain matrix
cornerNode = BVP.msh.cornerNode; % Corner nodes for monitoring
dof = 2; % Degrees of freedom per node
nn = 4; % Number of nodes per element

%% Initialize Arrays for Stresses and Strains
uniqueNodes = unique(et); % Unique nodes for averaging
cauchyStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % Cauchy stress
strainStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % Strain
pk1StressStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % PK1 stress

%% Loop Through Time Steps
for st = 1:size(u_steps, 2)
  cauchy = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3);
  pk1Stress = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3);
  strain = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3);
  u = u_steps(:, st); % Displacement for current time step
  
  % Loop Through Elements
  for e = 1:numberElements
    indice = et(e, :); % Node indices for the element
    elementDof = ed(e, :); % DOFs for the element
    elDisp = u(elementDof); % Element displacement
    elDisp = reshape(elDisp, dof, nn);
    c = 1; % Counter for Gauss points
    
    for j = 1:size(xik2, 2)
      for i = 1:size(xik1, 2)
        % Compute shape function derivatives in natural coordinates
        N1_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 1);
        N2_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 2);
        naturalDerivatives = [N1_d; N2_d]';
        
        % Compute Jacobian and derivatives in physical coordinates
        [~, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
        
        % Compute deformation gradient
        F = elDisp * XYDerivatives + eye(dof);
        
        % Compute strain tensor (Green-Lagrange)
        CC = F' * F;
        EE_tensor = 0.5 * (CC - eye(2));
        EE = voigt(EE_tensor); % Convert to Voigt notation
        
        % Compute stresses
        stress_2pk = C * E * voigt(EE_tensor); stress_2pk(3) = 2 * stress_2pk(3);
        stress_2pk_tensor = [stress_2pk(1) stress_2pk(3); stress_2pk(3) stress_2pk(2)];
        cauchy_tensor = F * stress_2pk_tensor * F';
        stress_1pk_tensor = inv(F) * cauchy_tensor;
        
        % Convert stresses to Voigt notation
        ccy = voigt(cauchy_tensor);
        pk1 = voigt(stress_1pk_tensor);
        
        % Store stresses and strains
        cauchy(e, c, :) = ccy;
        pk1Stress(e, c, :) = pk1;
        strain(e, c, :) = EE;
        c = c + 1;
      end
    end
    
  end
  
  % Average stresses and strains at nodes
  cauchy_xx_tmp = cauchy(:, :, 1); cauchy_xx = zeros(size(uniqueNodes, 1), 1);
  pk1Stress_xx_tmp = pk1Stress(:, :, 1); pk1Stress_xx = zeros(size(uniqueNodes, 1), 1);
  strain_xx_tmp = strain(:, :, 1); strain_xx = zeros(size(uniqueNodes, 1), 1);
  
  for ii = 1:size(uniqueNodes, 1)
    cauchy_xx(ii) = mean(cauchy_xx_tmp(ismember(et, uniqueNodes(ii))));
    pk1Stress_xx(ii) = mean(pk1Stress_xx_tmp(ismember(et, uniqueNodes(ii))));
    strain_xx(ii) = mean(strain_xx_tmp(ismember(et, uniqueNodes(ii))));
  end
  
  cauchyStep_xx(:, st) = cauchy_xx;
  pk1StressStep_xx(:, st) = pk1Stress_xx;
  strainStep_xx(:, st) = strain_xx;
end

%% Post-Processing
cornerStrechRatio = strainStep_xx(cornerNode, :);
cornerpk1Stress = pk1StressStep_xx(cornerNode, :);

cauchy_xx = cauchyStep_xx(:, end);
cauchy_xx_Surf = makeSurf(et, cauchy_xx);
pk1Stress_xx = pk1StressStep_xx(:, end);
pk1Stress_xx_Surf = makeSurf(et, pk1Stress_xx);

%% Assign Back to BVP
BVP.postproc.ST.cauchyStep_xx = cauchyStep_xx;
BVP.postproc.ST.pk1StressStep_xx = pk1StressStep_xx;
BVP.postproc.ST.strainStep_xx = strainStep_xx;
BVP.postproc.ST.cornerStrechRatio = cornerStrechRatio;
BVP.postproc.ST.cornerpk1Stress = cornerpk1Stress;
BVP.postproc.ST.cauchy_xx_Surf = cauchy_xx_Surf;
BVP.postproc.ST.cauchy_xx = cauchy_xx;
BVP.postproc.ST.pk1Stress_xx = pk1Stress_xx;
BVP.postproc.ST.pk1Stress_xx_Surf = pk1Stress_xx_Surf;

end
