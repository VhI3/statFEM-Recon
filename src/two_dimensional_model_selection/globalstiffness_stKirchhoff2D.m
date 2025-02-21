function [stiffness, R] = globalstiffness_stKirchhoff2D(BVP_tmp, u)
% Description:
% This function computes the global stiffness matrix and residual vector
% for a 2D finite element model using the St. Venant-Kirchhoff material model.
% It is designed to handle geometrically nonlinear elasticity problems.
%
% The function:
%   - Assembles the global stiffness matrix (including both material and geometric stiffness components)
%   - Computes the internal force (residual) vector
%   - Employs numerical integration with Gauss quadrature
%   - Utilizes the deformation gradient and Green-Lagrange strain tensor
%
% Inputs:
%   BVP_tmp - A structure containing:
%       • GDOFs            : Total global degrees of freedom
%       • nElm             : Number of elements
%       • nodeCoordinates  : Nodal coordinates
%       • elementNodes     : Element connectivity matrix
%       • elementDOFs      : Degrees of freedom per element
%       • GP               : Gauss points (xi1, xi2) and weights (w11, w12)
%       • C                : Constitutive matrix (plane strain/stress)
%       • E                : Young’s modulus
%
%   u       - Nodal displacement vector
%
% Outputs:
%   stiffness - Global stiffness matrix (GDOFs x GDOFs)
%   R         - Global residual vector (GDOFs x 1)
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon


%% Assign from BVP
GDOFs = BVP_tmp.GDOFs;
dof = 2; % degree of freedom
nElm = BVP_tmp.nElm;
nn = 4; % number of nodes for each element
nodeCoordinates = BVP_tmp.nodeCoordinates;
et = BVP_tmp.elementNodes;
ed = BVP_tmp.elementDOFs;
xi1 = BVP_tmp.GP.xi1;
xi2 = BVP_tmp.GP.xi2;
w1 = BVP_tmp.GP.w11;
w2 = BVP_tmp.GP.w12;
nnElem = 8;
C = BVP_tmp.C;
E = BVP_tmp.E;
%% Calculation
% initialazion
stiffness = zeros(GDOFs, GDOFs); % reserve stiffness matrix
R = zeros(GDOFs, 1); % reserve residual matrix
%
for e = 1:nElm % loop over elements
    indice = et(e, :); elementDof = ed(e, :); %  elementDof = 1:8;
    elDisp = u(elementDof);
    elDisp = reshape(elDisp, dof, nn);
    ke = zeros(size(elementDof, 2));
    re = zeros(size(elementDof, 2), 1);
    
    for j = 1:size(xi2, 2)
        
        for i = 1:size(xi1, 2)
            N1_d = Lagrange2D_d(xi1(i), xi2(j), 1, 1, [-1 1], [-1 1], 1);
            N2_d = Lagrange2D_d(xi1(i), xi2(j), 1, 1, [-1 1], [-1 1], 2);
            N1_d = N1_d([1 2 4 3]);
            N2_d = N2_d([1 2 4 3]);
            naturalDerivatives = [N1_d; N2_d]';
            %
            [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
            F = elDisp * XYDerivatives + eye(dof);
            CC = F' * F;
            EE_tensor = 0.5 * (CC - eye(2));
            stress_2pk = C * E * voigt(EE_tensor); stress_2pk(3) = 2 * stress_2pk(3);
            BN = zeros(3, nnElem);
            BG = zeros(4, nnElem);
            
            for k = 1:nn
              BN(:, k * 2 - 1:k * 2) = [F(1, 1) * XYDerivatives(k, 1) F(2, 1) * XYDerivatives(k, 1);
                F(1, 2) * XYDerivatives(k, 2) F(2, 2) * XYDerivatives(k, 2);
                F(1, 1) * XYDerivatives(k, 2) + F(1, 2) * XYDerivatives(k, 2) F(2, 1) * XYDerivatives(k, 2) + F(2, 2) * XYDerivatives(k, 1)];
              
              BG(:, k * 2 - 1:k * 2) = [XYDerivatives(k, 1) 0;
                XYDerivatives(k, 2) 0;
                0 XYDerivatives(k, 1);
                0 XYDerivatives(k, 2); ];
              
            end
            
            sigma = [stress_2pk(1) stress_2pk(3);
              stress_2pk(3) stress_2pk(2)];
            stan = zeros(4);
            stan(1:2, 1:2) = sigma;
            stan(3:4, 3:4) = sigma;
            %
            ke = ke + w1(i) * w2(j) * (BN' * C * E * BN + BG' * stan * BG) * det(JacobianMatrix);
            re = re + w1(i) * w2(j) * BN' * stress_2pk * det(JacobianMatrix);
        end
        
    end
    
    stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
    R(elementDof) = R(elementDof) + re;
end

end
