function [JacobianVector, invJacobian, XDerivatives] = Jacobian1D(nodeCoordinates, naturalDerivatives)
% JACOBIAN1D Computes the Jacobian matrix, its inverse, and global derivatives for 1D elements.
%
%   [JacobianMatrix, invJacobian, XDerivatives] = Jacobian1D(nodeCoordinates, naturalDerivatives)
%   calculates the Jacobian matrix for a 1D finite element, its inverse,
%   and the derivatives of the shape functions with respect to the global
%   coordinate (x).
%
% Inputs:
%   nodeCoordinates     - (nNodes x 1) Column vector of nodal coordinates [x].
%   naturalDerivatives  - (nNodes x 1) Column vector of derivatives of shape functions
%                         with respect to the natural coordinate ξ [dN/dξ].
%
% Outputs:
%   JacobianMatrix      - (1 x 1) Scalar Jacobian matrix (dx/dξ).
%   invJacobian         - (1 x 1) Scalar inverse of the Jacobian matrix (dξ/dx).
%   XDerivatives        - (nNodes x 1) Derivatives of shape functions with
%                         respect to global coordinates [dN/dx].
%
% Notes:
%   - This function is commonly used in finite element analysis to map natural
%     coordinates (ξ) to global coordinates (x).
%   - The Jacobian should not be close to zero to avoid numerical issues.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Step 1: Calculate the Jacobian matrix
% The Jacobian matrix maps natural coordinates (ξ) to global coordinates (x).
% J = dx/dξ
JacobianVector = naturalDerivatives*nodeCoordinates; % (1x1 scalar)
%% Step 2: Check for singularity
% Ensure the Jacobian is not too small to avoid numerical issues.
if abs(JacobianVector) < 1e-10
    error('Jacobian is close to singular. Check the element or mesh quality.');
end

%% Step 3: Compute the inverse of the Jacobian matrix
% Since this is 1D, the inverse is simply 1/J.
invJacobian = 1 / JacobianVector; % (1x1 scalar)

%% Step 4: Compute derivatives with respect to global coordinates
% Global derivatives: dN/dx = (dN/dξ) * (dξ/dx)
XDerivatives = naturalDerivatives * invJacobian; % (nNodes x 1)

end

