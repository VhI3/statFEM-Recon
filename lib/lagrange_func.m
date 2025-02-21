function [N] = lagrange_func(xi, xik, p)
% LAGRANGE_FUNC Computes Lagrange shape functions.
%
%   [N] = LAGRANGE_FUNC(xi, xik, p) evaluates the Lagrange shape functions
%   at specified points for a given polynomial degree and range.
%
% Inputs:
%   xi  - Points where the shape functions are evaluated (1D array).
%   xik - Nodes in the reference interval (1D array, typically [-1, 1]).
%   p   - Degree of the Lagrange polynomial.
%
% Outputs:
%   N   - Matrix of shape functions. N(i, j) is the value of the i-th shape
%         function at the j-th position in xi.
%
% Notes:
%   - This function implements the Lagrange shape functions using
%     the standard interpolation formula.
%   - Ensure xik is compatible with the degree p (length(xik) = p + 1).
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Last Modified: 18 February 2025

%% Input Validation
% Check if xik has the correct size for the given polynomial degree
if size(xik, 2) ~= (p + 1)
    error('The size of xik must match the polynomial degree p+1.');
end

%% Initialize Shape Functions
% Initialize the shape functions matrix
N = zeros(p + 1, size(xi, 2)); % (p + 1 shape functions for size(xi, 2) points)

%% Compute Shape Functions
% Loop over all xi values
for j = 1:size(xi, 2)
    % Loop over all shape functions
    for i = 1:(p + 1)
        N(i, j) = 1; % Initialize product for shape function i at xi(j)
        % Compute product for the Lagrange polynomial formula
        for k = 1:(p + 1)
            
            if k ~= i % Skip if k == i
              N(i, j) = N(i, j) * (xi(j) - xik(k)) / (xik(i) - xik(k));
            end
            
        end
        
    end
    
end

end
