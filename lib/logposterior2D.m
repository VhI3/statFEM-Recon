function log_posterior = logposterior2D(params, Cu, P, C_e, y_i, U_mean, x, n_sen, ii)
% LOGPOSTERIOR2D Computes the log-posterior for a 2D problem.
%   This function evaluates the log-posterior probability for a given set
%   of parameters, sensor data, and covariance matrices in a 2D analysis.
%
% Inputs:
%   params   - Hyperparameter vector [rho, length_scale, variance].
%              rho: Scaling factor.
%              length_scale: Length scale for the squared exponential kernel.
%              variance: Variance for the kernel.
%   Cu       - Covariance matrix for prior displacement.
%   P        - Projection matrix.
%   C_e      - Error covariance matrix.
%   y_i      - Sensor observations (n_sen x n_rep).
%   U_mean   - Mean displacement field.
%   x        - Sensor locations.
%   n_sen    - Number of sensors.
%   ii       - Index for specific functionality (e.g., selecting components).
%
% Outputs:
%   log_posterior - Computed log-posterior probability value.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Extract Parameters
rho = params(1); % Scaling factor (rho)

%% Initialize Variables
n_rep = size(y_i, 2); % Number of realizations (repetitions)

% Compute the covariance matrix
Sigma = rho^2 * P * Cu * P' + C_e + sqexp(x, x, params(2), params(3));

% Perform Cholesky decomposition
L = chol(Sigma, 'lower'); % Lower triangular Cholesky factor
invSigma = chol_solve(L, eye(n_sen)); % Efficient inversion using Cholesky factorization

%% Compute Log-Posterior Components
% Constant term
alpha1 = n_rep * n_sen * log(2 * pi);

% Log determinant term
alpha2 = 2 * n_rep * sum(log(diag(L)));

% Data fidelity term
alpha3 = 0;
for i = 1:n_rep
    residual = y_i(ii:2:end, i) - rho * P * U_mean; % Residual for current realization
    alpha3 = alpha3 + residual' * invSigma * residual; % Quadratic form
end

%% Compute Log-Posterior
log_posterior = 0.5 * (alpha1 + alpha2 + alpha3);

end
