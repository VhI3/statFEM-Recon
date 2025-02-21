function BVP = infinitePlate_statFEMbalancedStatus_ST(BVP)
% INFINITEPLATE_STATFEMBALANCEDSTATUS_ST
% Calculate the posterior residual forces for the St. Venant Kirchhoff model.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing preprocessed data.
%
% Outputs:
%   BVP - Updated structure with balanced equation results.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
dbc = BVP.BC.dbc;
prescribedDOFs = BVP.BC.prescribedDOFs;
force = BVP.BC.force;

% Temporary structure for element-specific computations
BVP_tmp.GP = BVP.proc.GP;
BVP_tmp.nElm = BVP.msh.nElm;
BVP_tmp.nodeCoordinates = BVP.msh.nodeCoordinates;
BVP_tmp.elementNodes = BVP.msh.elementNodes;
BVP_tmp.elementDOFs = BVP.msh.elementDOFs;
BVP_tmp.GDOFs = BVP.msh.GDOFs;

% Material properties
BVP_tmp.E = BVP.UQ.mu_E;
BVP_tmp.nu = BVP.material.nu;
BVP_tmp.lamda = BVP.material.lamda(BVP_tmp.E);
BVP_tmp.mu = BVP.material.mu(BVP_tmp.E);

% Displacement from StatFEM
u = BVP.statFEM.ST.mean_u_y;

%% Initialize Residual Force Vector
GDOFs = BVP_tmp.GDOFs;
R = zeros(GDOFs, 1);

%% Element-Level Residual Force Calculation
for e = 1:BVP_tmp.nElm
    % Element-specific data
    indice = BVP_tmp.elementNodes(e, :);
    elementDof = BVP_tmp.elementDOFs(e, :);
    elDisp = reshape(u(elementDof), 2, []); % Reshape displacement into 2 DOFs
    
    % Element residual forces
    re = zeros(length(elementDof), 1);
    
    % Gauss Quadrature Integration
    for j = 1:length(BVP_tmp.GP.xi2)
      
      for i = 1:length(BVP_tmp.GP.xi1)
        % Shape Function Derivatives
        N1_d = Lagrange2D_d(BVP_tmp.GP.xi1(i), BVP_tmp.GP.xi2(j), 1, 1, [-1 1], [-1 1], 1);
        N2_d = Lagrange2D_d(BVP_tmp.GP.xi1(i), BVP_tmp.GP.xi2(j), 1, 1, [-1 1], [-1 1], 2);
        naturalDerivatives = [N1_d([1 2 4 3]); N2_d([1 2 4 3])]';
        
        % Jacobian and Derivatives
        [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(BVP_tmp.nodeCoordinates(indice, :), naturalDerivatives);
        
        % Deformation Gradient (F)
        F = elDisp * XYDerivatives + eye(2);
        
        % Green-Lagrange Strain Tensor (E)
        CC = F' * F;
        EE = 0.5 * (CC - eye(2));
        
        % Stress Tensor (PK2)
        stress = BVP_tmp.lamda * trace(EE) * eye(2) + 2 * BVP_tmp.mu * EE;
        stressVec = [stress(1, 1); stress(2, 2); stress(1, 2)];
        
        % B-Matrix
        BN = zeros(3, 8);
        
        for k = 1:4
          BN(:, k * 2 - 1:k * 2) = [
            F(1, 1) * XYDerivatives(k, 1), F(2, 1) * XYDerivatives(k, 1);
            F(1, 2) * XYDerivatives(k, 2), F(2, 2) * XYDerivatives(k, 2);
            F(1, 1) * XYDerivatives(k, 2) + F(1, 2) * XYDerivatives(k, 1), ...
            F(2, 1) * XYDerivatives(k, 2) + F(2, 2) * XYDerivatives(k, 1)
            ];
        end
        
        % Residual Contribution
        re = re + BVP_tmp.GP.w11(i) * BVP_tmp.GP.w12(j) * BN' * stressVec * det(JacobianMatrix);
      end
      
    end
    
    % Assemble Element Residuals
    R(elementDof) = R(elementDof) + re;
end

%% Compute Balanced Forces
b = force - R;
b(prescribedDOFs) = dbc - u(prescribedDOFs);

% Separate Balanced Equations into x and y components
balancedEq_x = b(1:2:end);
balancedEq_y = b(2:2:end);

% Surface Representation
balancedEq_x_Surf = makeSurf(BVP_tmp.elementNodes, balancedEq_x);
balancedEq_y_Surf = makeSurf(BVP_tmp.elementNodes, balancedEq_y);

% Norm of Balanced Equation
norm_balancedEq = norm(b);

%% Assign Back to BVP
BVP.statFEM.ST.balancedEq = b;
BVP.statFEM.ST.balancedEq_x = balancedEq_x;
BVP.statFEM.ST.balancedEq_y = balancedEq_y;
BVP.statFEM.ST.balancedEq_x_Surf = balancedEq_x_Surf;
BVP.statFEM.ST.balancedEq_y_Surf = balancedEq_y_Surf;
BVP.statFEM.ST.norm_balancedEq = norm_balancedEq;

end
