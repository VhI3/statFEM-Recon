function BVP = infinitePlate_statFEM_postprocessST(BVP)
% INFINITEPLATE_STATFEM_POSTPROCESSST Computes stress distribution for St. Venant material in statFEM.
%
% This function calculates stresses at integration points and averages them
% to nodes, providing stress distributions and strain data for visualization
% and analysis in a St. Venant material model.
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
xik1 = [-1 1]; % Integration points in xi1 direction
xik2 = [-1 1]; % Integration points in xi2 direction
numberElements = BVP.msh.nElm; % Number of elements
u_steps = BVP.statFEM.ST.mean_u_y; % Posterior mean displacements
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
et = BVP.msh.elementNodes; % Element connectivity
ed = BVP.msh.elementDOFs; % Element degrees of freedom
E = BVP.UQ.mu_E; % Mean Young's modulus
C = BVP.material.C_pstrain; % Material stiffness matrix
dof = 2; % Degrees of freedom per node
nn = 4; % Number of nodes per element
cornerNode = BVP.msh.cornerNode; % Corner node indices

%% Initialize storage for stresses and strains
uniqueNodes = unique(et); % Unique node indices for averaging
cauchyStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % Storage for Cauchy stress
strainStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % Storage for strain
pk1StressStep_xx = zeros(size(uniqueNodes, 1), size(u_steps, 2)); % Storage for 1st PK stress

%% Loop over displacement steps
for st = 1:size(u_steps, 2)
    cauchy = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % Cauchy stress tensor
    pk1Stress = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % 1st PK stress tensor
    strain = zeros(numberElements, size(xik1, 2) * size(xik2, 2), 3); % Strain tensor
    u = u_steps(:, st); % Displacement for current step
    
    % Loop over elements
    for e = 1:numberElements
        indice = et(e, :); % Nodes of the element
        elementDof = ed(e, :); % DOFs of the element
        elDisp = u(elementDof); % Element displacements
        elDisp = reshape(elDisp, dof, nn); % Reshape to 2x4 matrix
        c = 1; % Counter for integration points
        
        % Loop over integration points
        for j = 1:size(xik2, 2)
          
          for i = 1:size(xik1, 2)
            % Compute derivatives of shape functions in natural coordinates
            N1_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 1);
            N2_d = Lagrange2D_d(xik1(i), xik2(j), 1, 1, xik1, xik2, 2);
            naturalDerivatives = [N1_d; N2_d]';
            
            % Compute Jacobian matrix and derivatives in global coordinates
            [~, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
            
            % Compute deformation gradient (F)
            F = elDisp * XYDerivatives + eye(dof);
            CC = F' * F; % Right Cauchy-Green tensor
            EE_tensor = 0.5 * (CC - eye(dof)); % Green-Lagrange strain tensor
            EE = voigt(EE_tensor); % Strain in Voigt notation
            
            % Compute 2nd Piola-Kirchhoff stress
            stress_2pk = C * E * EE;
            stress_2pk(3) = 2 * stress_2pk(3); % Correct shear stress
            stress_2pk_tensor = [stress_2pk(1) stress_2pk(3); stress_2pk(3) stress_2pk(2)];
            
            % Compute Cauchy stress
            cauchy_tensor = F * stress_2pk_tensor * F';
            ccy = voigt(cauchy_tensor); % Cauchy stress in Voigt notation
            
            % Compute 1st Piola-Kirchhoff stress
            stress_1pk_tensor = inv(F) * cauchy_tensor;
            pk1 = voigt(stress_1pk_tensor); % 1st PK stress in Voigt notation
            
            % Store results
            cauchy(e, c, :) = ccy;
            pk1Stress(e, c, :) = pk1;
            strain(e, c, :) = EE;
            c = c + 1; % Increment counter
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
    
    % Store step results
    cauchyStep_xx(:, st) = cauchy_xx;
    pk1StressStep_xx(:, st) = pk1Stress_xx;
    strainStep_xx(:, st) = strain_xx;
end

% Compute corner stretch ratios and stresses
cornerStrechRatio = strainStep_xx(cornerNode, :) + 1;
cornerpk1Stress = pk1StressStep_xx(cornerNode, :);

% Final step visualization
cauchy_xx = cauchyStep_xx(:, end);
cauchy_xx_Surf = makeSurf(et, cauchy_xx);
pk1Stress_xx = pk1StressStep_xx(:, end);
pk1Stress_xx_Surf = makeSurf(et, pk1Stress_xx);

%% Assign results back to BVP
BVP.statFEM.postproc.ST.cauchyStep_xx = cauchyStep_xx;
BVP.statFEM.postproc.ST.pk1StressStep_xx = pk1StressStep_xx;
BVP.statFEM.postproc.ST.strainStep_xx = strainStep_xx;
BVP.statFEM.postproc.ST.cornerStrechRatio = cornerStrechRatio;
BVP.statFEM.postproc.ST.cornerpk1Stress = cornerpk1Stress;
BVP.statFEM.postproc.ST.cauchy_xx_Surf = cauchy_xx_Surf;
BVP.statFEM.postproc.ST.cauchy_xx = cauchy_xx;
BVP.statFEM.postproc.ST.pk1Stress_xx = pk1Stress_xx;
BVP.statFEM.postproc.ST.pk1Stress_xx_Surf = pk1Stress_xx_Surf;

end
