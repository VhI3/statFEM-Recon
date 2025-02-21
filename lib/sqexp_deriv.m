function deriv = sqexp_deriv(x1, x2, sigma, l)
% SQEXP_DERIV Computes the derivatives of the squared exponential covariance matrix.
%
%   deriv = SQEXP_DERIV(x1, x2, sigma, l) computes the derivatives of the
%   squared exponential covariance matrix with respect to sigma and l.
%
% Inputs:
%   x1    - (n1 x d) Matrix of input points, where each row is a d-dimensional point.
%   x2    - (n2 x d) Matrix of input points, where each row is a d-dimensional point.
%   sigma - Logarithm of the signal variance (scalar).
%   l     - Logarithm of the length scale (scalar).
%
% Outputs:
%   deriv - (n1 x n2 x 2) Tensor containing derivatives:
%           deriv(:,:,1) = d/d(sigma) of the squared exponential covariance.
%           deriv(:,:,2) = d/d(l) of the squared exponential covariance.
%
% Notes:
%   - This function calculates the derivatives of the squared exponential kernel:
%     k(x, y) = exp(2*sigma) * exp(-0.5 * ||x - y||^2 * exp(-2*l)).
%   - These derivatives are useful for hyperparameter optimization in Gaussian Process Regression.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Adapted from the statFEM Python library.

%% Compute Pairwise Squared Distances
r2 = pdist2(x1, x2) .^ 2; % Pairwise squared Euclidean distances

%% Compute Covariance Matrix
K = sqexp(x1, x2, sigma, l); % Compute squared exponential covariance

%% Compute Derivatives
deriv = zeros(size(x1,1), size(x2,1), 2);
deriv(:,:,1) = 2 * K;               % Derivative with respect to sigma
deriv(:,:,2) = K .* r2 .* exp(-2*l); % Derivative with respect to l
end
