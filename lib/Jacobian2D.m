function [JacobianMatrix, invJacobian, XYDerivatives] = Jacobian2D(nodeCoordinates, naturalDerivatives)
% JACOBIAN2D Computes the Jacobian matrix, its inverse, and global derivatives.
%
%   [JacobianMatrix, invJacobian, XYDerivatives] = Jacobian2D(nodeCoordinates, naturalDerivatives)
%   calculates the Jacobian matrix for a 2D finite element, its inverse,
%   and the derivatives of the shape functions with respect to global
%   coordinates (x, y).
%
% Inputs:
%   nodeCoordinates     - (nNodes x 2) Matrix of nodal coordinates [x, y].
%   naturalDerivatives  - (nNodes x 2) Matrix of derivatives of shape functions
%                         with respect to natural coordinates [dN/dξ, dN/dη].
%
% Outputs:
%   JacobianMatrix      - (2 x 2) Jacobian matrix.
%   invJacobian         - (2 x 2) Inverse of the Jacobian matrix.
%   XYDerivatives       - (nNodes x 2) Derivatives of shape functions with
%                         respect to global coordinates [dN/dx, dN/dy].
%
% Notes:
%   - This function is commonly used in finite element analysis to map natural
%     coordinates (ξ, η) to global coordinates (x, y).
%   - The determinant of the Jacobian matrix should not be close to zero
%     to avoid singularities.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Step 1: Calculate the Jacobian matrix
% The Jacobian matrix maps natural to global coordinates.
% J = [∂x/∂ξ ∂x/∂η; ∂y/∂ξ ∂y/∂η]
JacobianMatrix = nodeCoordinates' * naturalDerivatives;

%% Step 2: Check for singularity
% Ensure the determinant of the Jacobian is not close to zero
if abs(det(JacobianMatrix)) < 1e-10
    error('Jacobian is close to singular. Check the element or mesh quality.');
end

%% Step 3: Calculate the inverse of the Jacobian matrix
% Using backslash operator (\) for numerical stability instead of inv().
invJacobian = JacobianMatrix \ eye(2);

%% Step 4: Calculate derivatives with respect to global coordinates
% Global derivatives: [dN/dx dN/dy] = [dN/dξ dN/dη] * (Jacobian inverse)
XYDerivatives = naturalDerivatives * invJacobian;
end
