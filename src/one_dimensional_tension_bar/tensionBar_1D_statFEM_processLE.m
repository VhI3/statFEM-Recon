function BVP = tensionBar_1D_statFEM_processLE(BVP)
% This function performs the statistical finite element method (statFEM)
% analysis for a 1D tension bar using linear elasticity. The analysis
% incorporates uncertainty quantification by updating the displacement
% field based on observed noisy mismatched data and prior covariance structures.
%
% The function computes:
%   - Posterior mean displacement vector (mu_u_y)
%   - Posterior covariance matrix (C_u_y)
%   - Credible intervals for displacements (ci_u_y)
%   - Posterior mean and covariance of true response (mu_z, C_z)
%
% Method:
% The posterior estimates are computed using Bayesian inference with
% Gaussian priors and likelihoods. Cholesky decomposition is employed
% for numerical stability and efficient matrix inversion.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing:
%       • Material properties: sig_E
%       • Prior mean and covariance: mu_u_active, C_u_active
%       • Projection matrix: P_active
%       • Observation covariance: C_e
%       • Observation data: y_obs
%       • Model mismatch covariance: C_d
%       • Hyperparameters
%
% Outputs:
%   BVP - Updated structure containing:
%       • Posterior mean displacement (mu_u_y)
%       • Posterior covariance (C_u_y)
%       • 95% credible intervals (ci_u_y)
%       • Posterior mean measurement predictions (mu_z)
%       • Posterior covariance of measurements (C_z)
%       • 95% credible intervals of true response (ci_z)
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Assign the parameters
GDof = BVP.fem.GDof;
sig_E = BVP.material.sig_E;
C_u_active = BVP.UQ.LE.C_u_active_pc;
mu_u_active = BVP.UQ.LE.mu_u_active_pc;
activeDOFs = BVP.fem.activeDOFs;
C_d = BVP.statFEM.C_d;
rho = BVP.statFEM.rho;
C_e = BVP.obs.C_e;
P_active = BVP.obs.P_active;
y_obs = BVP.obs.y_obs;
nrep = BVP.obs.nrep;
%% Initialize Posterior Mean and Covariance
mu_u_y = zeros(GDof, 1);
C_u_y = zeros(GDof, GDof);

% Stabilize covariance matrix
C_u_active = C_u_active +1e-1 * sig_E * eye(size(activeDOFs, 2));

%% Compute (C_d + C_e) Inverse using Cholesky Decomposition
L_CdCe = chol(C_d + C_e, 'lower'); % Cholesky factorization: C_d + C_e = L_CdCe * L_CdCe'
inv_CdCe = chol_solve(L_CdCe, eye(size(C_d + C_e))); %#ok<ELARLOG> % Efficient inversion: inv(C_d + C_e)

%% Compute Inverse of C_u_active using Cholesky Decomposition
L_Cu = chol(C_u_active, 'lower'); % C_u_active = L_Cu * L_Cu'
inv_Cu = chol_solve(L_Cu, eye(size(C_u_active))); % inv(C_u_active)

%% Compute the Posterior Covariance (C_u_y_active)
% Original: C_u_y_active = inv(rho^2 * nrep * P_active' * inv(C_d + C_e) * P_active + inv(C_u_active));
M = rho ^ 2 * nrep * (P_active' * inv_CdCe * P_active) + inv_Cu;
L_M = chol(M, 'lower'); % Cholesky factor of M
C_u_y_active = chol_solve(L_M, eye(size(M))); % Inverse using Cholesky

%% Compute the Posterior Mean (mu_u_y_active)
% Original: mu_u_y_active = C_u_y_active * (rho * P_active' * inv(C_d + C_e) * sum(y_obs, 2) + inv(C_u_active) * mu_u_active);

right_hand_side = rho * (P_active' * inv_CdCe * sum(y_obs(:, 1:nrep), 2)) + (inv_Cu * mu_u_active);
mu_u_y_active = C_u_y_active * right_hand_side;

% Mean vector and Covariance matrix of posterior
mu_u_y(2:end) = rho * mu_u_y_active;
C_u_y(2:end, 2:end) = rho ^ 2 * C_u_y_active;
ci_u_y = sqrt(diag(C_u_y)) * 1.96;

% True response
mu_z = rho * P_active * mu_u_y_active;
C_z = rho ^ 2 * P_active * C_u_y_active * P_active' + C_d;
ci_z = sqrt(diag(C_z)) * 1.96;

%% Assign back to BVP
BVP.statFEM.mu_u_y = mu_u_y;
BVP.statFEM.C_u_y = C_u_y;
BVP.statFEM.ci_u_y = ci_u_y;
BVP.statFEM.mu_z = mu_z;
BVP.statFEM.C_z = C_z;
BVP.statFEM.ci_z = ci_z;
end
