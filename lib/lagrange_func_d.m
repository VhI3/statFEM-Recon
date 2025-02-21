function [Nd] = lagrange_func_d(xi, xik, p)
% LAGRANGE_FUNC_D Computes the derivative of Lagrange shape functions.
%
%   [Nd] = LAGRANGE_FUNC_D(xi, xik, p) evaluates the derivatives of the
%   Lagrange shape functions at specified points for a given polynomial degree.
%
% Inputs:
%   xi  - Points where the shape function derivatives are evaluated (1D array).
%   xik - Nodes in the reference interval (1D array, typically [-1, 1]).
%   p   - Degree of the Lagrange polynomial.
%
% Outputs:
%   Nd  - Matrix of derivative values. Nd(i, j) is the derivative of the i-th
%         shape function at the j-th position in xi.
%
% Notes:
%   - This function computes the derivatives using the standard Lagrange
%     polynomial derivative formula.
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

%% Initialize Derivative Matrix
% Initialize the matrix for storing the derivatives
Nd = zeros(p + 1, size(xi, 2)); % (p + 1 shape functions for size(xi, 2) points)

%% Compute Derivatives
% Loop over all xi values
for j = 1:size(xi, 2)
    % Loop over all shape functions
    for i = 1:(p + 1)
        Nd(i, j) = 0; % Initialize derivative value for the i-th shape function at xi(j)
        
        for m = 1:(p + 1)
            
            if m ~= i % Skip m == i
                temp = 1; % Initialize temporary product for the derivative formula
                % Compute the product excluding indices k == i and k == m
                for k = 1:(p + 1)
                    
                    if k ~= i && k ~= m % Skip k == i and k == m
                        temp = temp * (xi(j) - xik(k)) / (xik(k) - xik(i));
                    end
                    
                end
                
                % Accumulate the derivative contribution from the current term
                Nd(i, j) = Nd(i, j) - temp / (xik(m) - xik(i));
            end
            
        end
        
    end
    
end

end
