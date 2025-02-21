function cov = sqexp(x1, x2, sigma, l)
% SQEXP Computes the squared exponential covariance matrix.
%
%   cov = SQEXP(x1, x2, sigma, l) computes the squared exponential
%   covariance matrix for input points x1 and x2.
%
% Inputs:
%   x1    - (n1 x d) Matrix of input points, where each row is a d-dimensional point.
%   x2    - (n2 x d) Matrix of input points, where each row is a d-dimensional point.
%   sigma - Logarithm of the signal variance (scalar).
%   l     - Logarithm of the length scale (scalar).
%
% Outputs:
%   cov - (n1 x n2) Covariance matrix computed using the squared exponential kernel.
%
% Notes:
%   - The squared exponential kernel is defined as:
%     k(x, y) = exp(2*sigma) * exp(-0.5 * ||x - y||^2 * exp(-2*l)).
%   - This kernel is commonly used in Gaussian process regression.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Adapted from the statFEM Python library.

%% Compute Pairwise Squared Distances
r2 = pdist2(x1, x2) .^ 2; % Pairwise squared Euclidean distances

%% Compute Covariance Matrix
% The squared exponential kernel formula
cov = exp(2 * sigma) * exp(-0.5 * r2 * exp(-2 * l));
end
