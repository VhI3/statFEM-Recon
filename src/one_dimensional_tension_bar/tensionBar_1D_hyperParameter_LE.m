function BVP = tensionBar_1D_hyperParameter_LE(BVP)
% TENSIONBAR_1D_HYPERPARAMETER_LE - Estimates hyperparameters for StatFEM in 1D tension bar problem.
%
%   This function estimates the hyperparameters (ρ, σ_d, l_d) for the 1D
%   tension bar problem using a log-likelihood optimization approach.
%
%   The optimization is performed using `fminunc` with the Trust-Region
%   method, leveraging analytical gradients.
%
% Inputs:
%   BVP - A structure containing:
%     - BVP.obs.rho_sample      : Initial scaling factor ρ.
%     - BVP.obs.sig_d_sample    : Initial standard deviation σ_d.
%     - BVP.obs.l_d_sample      : Initial length scale l_d.
%     - BVP.UQ.LE.C_u_active_pc : Covariance matrix of the active PCE response.
%     - BVP.obs.C_e             : Measurement noise covariance matrix.
%     - BVP.obs.P_active        : Projection matrix for observations.
%     - BVP.UQ.LE.mu_u_active_pc: Mean displacement from PCE.
%     - BVP.obs.senCoor         : Sensor coordinates.
%     - BVP.obs.y_obs           : Experimental displacement observations.
%     - BVP.obs.nrep            : Number of experiment repetitions.
%     - BVP.obs.nsen            : Number of sensors.
%     - BVP.obs.S               : Non-linearity indicator.
%
% Outputs:
%   BVP - Updated structure with identified hyperparameters:
%     - BVP.statFEM.rho   : Estimated scaling factor ρ.
%     - BVP.statFEM.sig_d : Estimated standard deviation σ_d.
%     - BVP.statFEM.l_d   : Estimated length scale l_d.
%
% Optimization:
%   - Uses `fminunc` with the **trust-region** method.
%   - Employs **analytical gradients** for improved convergence.
%   - Finite difference gradients (central differences) for robustness.
%
% Example:
%   BVP = tensionBar_1D_hyperParameter_LE(BVP);
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% **Step 1: Assign Input Parameters from BVP**
rho_sample = BVP.obs.rho_sample;
sig_d_sample = BVP.obs.sig_d_sample;
l_d_sample = BVP.obs.l_d_sample;
C_u_active_pc = BVP.UQ.LE.C_u_active_pc;
C_e = BVP.obs.C_e;
P_active = BVP.obs.P_active;
mu_u_active_pc = BVP.UQ.LE.mu_u_active_pc;
senCoor = BVP.obs.senCoor;
y_obs = BVP.obs.y_obs;
nrep = BVP.obs.nrep;
nsen = BVP.obs.nsen;
S = BVP.obs.S; % Non-linearity indicator

%% **Step 2: Set Initial Guess for Optimization**
start_w = [1, log(0.8), log(5)]; % Initial parameter guesses

% Define negative log-likelihood function
negLogLikelihoodFunc = @(w) negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);

%% **Step 3: Optimization Options**
options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ... % Use trust-region algorithm
  'SpecifyObjectiveGradient', true, ... % Enable analytical gradients
  'FiniteDifferenceType', 'central', ... % Use central finite differences for numerical stability
  'Diagnostics', 'on', ... % Enable diagnostics
  'Display', 'iter', ... % Show iterative optimization progress
  'MaxIterations', 200, ... % Set iteration limit
  'OptimalityTolerance', 1e-6, ...
  'StepTolerance', 1e-6);

%% **Step 4: Gradient Checking (Optional)**
% Uncomment the following line to check gradients before optimization:
% checkGradients(negLogLikelihoodFunc, start_w);

%% **Step 5: Run the Optimizer**
[identified_w] = fminunc(negLogLikelihoodFunc, start_w, options);

%% **Step 6: Display Identified Hyperparameters**
if S == 0 % Reference check for linear case
  fprintf('\nActual Hyperparameters:\n');
  fprintf('ρ = %f,\t σ_d = %f,\t l_d = %f \n', rho_sample, sig_d_sample, l_d_sample);
end

fprintf('\nEstimated Hyperparameters:\n');
fprintf('ρ = %f,\t σ_d = %f,\t l_d = %f \n', identified_w(1), exp(identified_w(2)), exp(identified_w(3)));

C_d = sqexp(senCoor, senCoor, identified_w(2), identified_w(3));
%% **Step 7: Assign Results to BVP Structure**
BVP.statFEM.rho = identified_w(1);
BVP.statFEM.sig_d = exp(identified_w(2));
BVP.statFEM.l_d = exp(identified_w(3));
BVP.statFEM.C_d = C_d;
end

function [f, g] = negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
% NEGATIVELOGLIKELIHOOD - Computes the negative log-likelihood and its gradient.
%
%   This function evaluates the negative log-posterior (log-likelihood) for
%   hyperparameter optimization in the 1D tension bar problem.
%
% Inputs:
%   w               - Vector of hyperparameters [ρ, log(σ_d), log(l_d)].
%   C_u_active_pc   - Covariance matrix of the active PCE response.
%   P_active        - Projection matrix for observations.
%   C_e             - Measurement noise covariance matrix.
%   y_obs           - Experimental observations (m sensors, k repetitions).
%   mu_u_active_pc  - Mean displacement from PCE.
%   senCoor         - Sensor coordinates.
%   nrep            - Number of experiment repetitions.
%   nsen            - Number of sensors.
%
% Outputs:
%   f - Scalar value of the negative log-likelihood.
%   g - (3x1) Gradient vector of the negative log-likelihood with respect to w.
%
% Usage Example:
%   [f, g] = negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);
%
% Notes:
%   - The function is used in optimization (e.g., `fminunc`) to estimate optimal
%     hyperparameters for the StatFEM approach.
%   - Computes both function value and gradient to improve optimization efficiency.
%
% Author: Vahab Narouie, TU-Braunschweig, 2024

%% **Step 1: Compute Negative Log-Posterior (Objective Function)**
f = logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);

%% **Step 2: Compute Gradient (Jacobian) of Log-Posterior**
g = logpost_deriv1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);
end

function log_posterior = logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
% LOGPOSTERIOR1D - Computes the negative log-posterior for hyperparameter estimation.
%
%   This function calculates the log-posterior (negative log-likelihood)
%   for hyperparameter optimization in the 1D tension bar problem using a
%   squared exponential covariance model.
%
% Inputs:
%   w               - Vector of hyperparameters [ρ, log(σ_d), log(l_d)].
%   C_u_active_pc   - Covariance matrix of the active PCE response.
%   P_active        - Projection matrix for observations.
%   C_e             - Measurement noise covariance matrix.
%   y_obs           - Experimental observations (m sensors, k repetitions).
%   mu_u_active_pc  - Mean displacement from PCE.
%   senCoor         - Sensor coordinates.
%   nrep            - Number of experiment repetitions.
%   nsen            - Number of sensors.
%
% Outputs:
%   log_posterior - Scalar value of the negative log-posterior.
%
% Usage Example:
%   log_posterior = logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);
%
% Notes:
%   - Used in hyperparameter optimization (e.g., `fminunc`).
%   - Efficiently computes Cholesky decomposition for numerical stability.
%   - Implements Gaussian process modeling with squared exponential kernel.
%
% Author: Vahab Narouie, TU-Braunschweig, 2024

%% **Step 1: Extract Hyperparameters**
rho = w(1);
sig_d = w(2);
l_d = w(3);

%% **Step 2: Construct Covariance Matrix**
Sigma = rho ^ 2 * P_active * C_u_active_pc * P_active' + C_e + sqexp(senCoor, senCoor, sig_d, l_d);

%% **Step 3: Cholesky Decomposition for Numerical Stability**
L = chol(Sigma, 'lower'); % Lower triangular Cholesky factor
invSigma = chol_solve(L, eye(nsen)); % Efficient inversion using Cholesky factorization

%% **Step 4: Compute Log-Likelihood Components**
alpha1 = nrep * nsen * log(2 * pi); % Constant term
alpha2 = 2 * nrep * sum(log(diag(L))); % Log determinant term

% Data fidelity term
alpha3 = 0;

for i = 1:nrep
  residual = y_obs(:, i) - rho * P_active * mu_u_active_pc; % Compute residual
  alpha3 = alpha3 + residual' * invSigma * residual; % Quadratic form
end

%% **Step 5: Compute Negative Log-Posterior**
log_posterior = 0.5 * (alpha1 + alpha2 + alpha3);
end

function deriv = logpost_deriv1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
% LOGPOST_DERIV1D - Computes the derivative of the log-posterior function
%
% Computes the gradient of the log-posterior w.r.t hyperparameters:
%   1. Scaling factor (rho)
%   2. Standard deviation of model error (sigma_d)
%   3. Correlation length (l_d)
%
% Inputs:
%   w                   - [rho, log(sigma_d), log(l_d)] hyperparameters
%   C_u_active_pc       - Covariance matrix of active Displacement DOFs
%   P_active            - Projection matrix
%   C_e                 - Observation noise covariance
%   y_obs               - Experimental observations (n_sen × n_rep)
%   mu_u_active_pc      - Active Mean displacement from PC expansion
%   senCoor             - Sensor coordinates
%   nrep                - Number of realizations
%   nsen                - Number of sensors
%
% Output:
%   deriv  - Gradient vector of log-posterior w.r.t hyperparameters
%
% Author: Vahab Narouie, TU-Braunschweig, 2024

%% Step 1: Extract Hyperparameters
rho = w(1);
sig_d = w(2);
l_d = w(3);
%
Sigma = rho ^ 2 * P_active * C_u_active_pc * P_active' + C_e + sqexp(senCoor, senCoor, sig_d, l_d);
% Perform Cholesky decomposition
L = chol(Sigma, 'lower'); % Lower triangular Cholesky factor
invSigma = chol_solve(L, eye(nsen)); % Efficient inversion using Cholesky factorization

% Compute derivatives of the covariance function
K_deriv = sqexp_deriv(senCoor, senCoor, sig_d, l_d);
Cd_sigd = K_deriv(:, :, 1); % ∂C_d / ∂σ_d
Cd_ld = K_deriv(:, :, 2); % ∂C_d / ∂l_d

%% Step 3: Initialize Derivative Vector
deriv = zeros(3, 1);

%% Step 4: Compute Derivative w.r.t ρ (rho)
deriv(1) = nrep * rho * 2 * trace(invSigma * P_active * C_u_active_pc * P_active');

for i = 1:nrep
  residual = y_obs(:, i) - rho * P_active * mu_u_active_pc; % Compute residual
  a1 = (P_active * mu_u_active_pc)' * invSigma * residual;
  a2 = residual' * invSigma * (rho * 2 * P_active * C_u_active_pc * P_active' * invSigma) * residual;
  a3 = residual' * invSigma * P_active * mu_u_active_pc;
  deriv(1) = deriv(1) - (a1 + a2 + a3);
end

deriv(1) = 0.5 * deriv(1);

%% Step 5: Compute Derivative w.r.t σ_d
deriv(2) = nrep * trace(invSigma * Cd_sigd);

for i = 1:nrep
  residual = y_obs(:, i) - rho * P_active * mu_u_active_pc;
  deriv(2) = deriv(2) - residual' * invSigma * Cd_sigd * invSigma * residual;
end

deriv(2) = 0.5 * deriv(2);

%% Step 6: Compute Derivative w.r.t l_d
deriv(3) = nrep * trace(invSigma * Cd_ld);

for i = 1:nrep
  residual = y_obs(:, i) - rho * P_active * mu_u_active_pc;
  deriv(3) = deriv(3) - residual' * invSigma * Cd_ld * invSigma * residual;
end

deriv(3) = 0.5 * deriv(3);
end
