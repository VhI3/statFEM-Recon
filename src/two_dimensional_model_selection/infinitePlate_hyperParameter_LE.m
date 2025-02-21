function BVP = infinitePlate_hyperParameter_LE(BVP)
% INFINITEPLATE_HYPERPARAMETER_LE Estimates hyperparameters for linear elasticity.
%   This function estimates the hyperparameters (rho, sig_d, l_d) for each
%   dimension (x and y) separately and computes corresponding covariance matrices.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       Observation data, FEM results, and initial hyperparameters.
%
% Outputs:
%   BVP - Updated BVP structure with estimated hyperparameters and covariance matrices.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
% Observational data
sig_d_sample = BVP.obs.sig_d_sample; % Initial sig_d sample
l_d_sample = BVP.obs.l_d_sample; % Initial l_d sample
rho_sample = BVP.obs.rho_sample; % Initial rho sample
nSen = BVP.obs.nSen; % Number of sensors
epsi = BVP.obs.epsi; % Noise level
senCoor = BVP.obs.senCoor; % Sensor coordinates
Y_exp = BVP.obs.Y_exp; % Experimental observations
P_active = BVP.obs.P_active; % Observation matrices for active DOFs

% FEM results
mean_u_active = BVP.lsNIPCE.LE.mean_u_active_pc; % Mean active displacements
cov_u_active = BVP.lsNIPCE.LE.cov_u_active_pc; % Covariance of active displacements

% Mesh details
DOFs = BVP.msh.DOFs; % Degrees of freedom per node

%% Initialize Covariance and Parameters
C_e = epsi * eye(nSen); % Error covariance matrix
startPar = [rho_sample, log(1), log(1)]; % Initial parameter guesses

% Storage for hyperparameters and covariance matrices
rho_est = zeros(DOFs, 1);
sigd_est = zeros(DOFs, 1);
ld_est = zeros(DOFs, 1);
C_d = cell(DOFs, 1);

%% Optimization for Hyperparameter Estimation
for ii = 1:DOFs
  % Define the objective function for optimization
  fun = @(x) logFunc2D(x, cov_u_active, P_active{ii}, C_e, Y_exp, mean_u_active, senCoor, nSen, ii);
  
  % Optimize using fminunc
  params_est = fminunc(fun, startPar);
  
  % Extract hyperparameters for this dimension
  rho_est(ii) = params_est(1);
  sigd_est(ii) = exp(params_est(2));
  ld_est(ii) = exp(params_est(3));
  
  % Compute the covariance matrix for this dimension
  C_d{ii} = sqexp(senCoor, senCoor, params_est(2), params_est(3));
end

%% Display Actual and Estimated Parameters
fprintf('\nActual parameters: rho = %f, sig_d = %f, l_d = %f\n', rho_sample, sig_d_sample, l_d_sample);

for ii = 1:DOFs
  fprintf('Dimension %d Estimated parameters: rho = %f, sig_d = %f, l_d = %f\n', ii, rho_est(ii), sigd_est(ii), ld_est(ii));
end

%% Assign Results Back to BVP
BVP.hyp.LE.sigd = sigd_est; % Estimated sig_d for each dimension
BVP.hyp.LE.ld = ld_est; % Estimated l_d for each dimension
BVP.hyp.LE.rho = rho_est; % Estimated rho for each dimension
BVP.hyp.LE.C_e = C_e; % Error covariance matrix
BVP.hyp.LE.C_d = C_d; % Covariance matrices for each dimension

end
