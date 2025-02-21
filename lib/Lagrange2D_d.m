function [N] = Lagrange2D_d(xi1, xi2, p1, p2, xik1, xik2, dim)
% LAGRANGE2D_D Computes the derivative of 2D Lagrange shape functions.
%
%   [N] = LAGRANGE2D_D(xi1, xi2, p1, p2, xik1, xik2, dim) evaluates the
%   derivatives of the Lagrange shape functions for 2D elements based on
%   the polynomial degrees and ranges in two directions.
%
% Inputs:
%   xi1  - Point(s) in the first direction (X or X1).
%   xi2  - Point(s) in the second direction (Y or X2).
%   p1   - Polynomial degree in the first direction.
%   p2   - Polynomial degree in the second direction.
%   xik1 - Range of xi1 (nodes in the reference interval, e.g., [-1, 1]).
%   xik2 - Range of xi2 (nodes in the reference interval, e.g., [-1, 1]).
%   dim  - Dimension for differentiation:
%          1 = Differentiate with respect to xi1 (X or X1),
%          2 = Differentiate with respect to xi2 (Y or X2).
%
% Outputs:
%   N    - Derivative vector of shape functions at positions xi1 and xi2.
%
% Notes:
%   - This function combines the derivatives in one direction and shape
%     functions in the other direction for 2D problems.
%   - Supports higher polynomial degrees for Lagrange interpolation.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Last Modified: 18 February 2025

%% Initialize Derivative Matrix
% Total number of shape functions: (p1+1) * (p2+1)
N = zeros(1, (p1 + 1) * (p2 + 1));

%% Compute Shape Functions and Their Derivatives
if dim == 1
    % Differentiate with respect to xi1
    N1 = lagrange_func_d(xi1, xik1, p1); % Derivatives in xi1 direction
    N2 = lagrange_func(xi2, xik2, p2); % Shape functions in xi2 direction
else
    % Differentiate with respect to xi2
    N1 = lagrange_func(xi1, xik1, p1); % Shape functions in xi1 direction
    N2 = lagrange_func_d(xi2, xik2, p2); % Derivatives in xi2 direction
end

%% Compute Derivatives for Each Shape Function
c = 1; % Counter for shape functions

for j = 1:(p2 + 1) % Loop over second dimension
    
    for i = 1:(p1 + 1) % Loop over first dimension
        % Combine derivatives and shape functions
        N(1, c) = N1(i, 1) * N2(j, 1);
        c = c + 1;
    end
    
end

end

%% TODO:
% 1. Validate for higher degrees (p1 = 2, 3; p2 = 2, 3) with intervals [-1, 1].
% 2. Check with different ranges: [0, 1], [-1, 1] for xi1 and xi2.
% 3. Extend tests for mixed ranges (e.g., xi1 in [-1, 1], xi2 in [0, 1]).
