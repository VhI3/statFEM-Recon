function BVP = infinitePlate_BayesFactor(BVP)
% INFINITEPLATE_BAYESFACTOR Computes the Bayes Factor for model comparison.
%
% This function calculates the Bayes Factor (BF) to compare the
% St. Venant (ST) material model and the Linear Elasticity (LE) model
% using observed data and prior information, considering separate
% hyperparameters for X and Y directions.
%
% Inputs:
%   BVP - Boundary value problem structure containing observed data,
%         model predictions, and hyperparameters for LE and ST models.
%
% Outputs:
%   BVP - Updated BVP structure containing Bayes Factor and log-likelihood
%         functions for LE and ST models.
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Assign from BVP
nSen = BVP.obs.nSen; % Number of sensors
nRep = BVP.obs.nSam; % Number of data samples (repetitions)

P_active = BVP.obs.P_active; % Active projection matrices for each DOF
Y_exp = BVP.obs.Y_exp; % Observed data

% Linear Elasticity (LE) model parameters
C_e_LE           = BVP.hyp.LE.C_e;
C_d_LE_x         = BVP.hyp.LE.C_d{1}; % C_d for X-direction
C_d_LE_y         = BVP.hyp.LE.C_d{2}; % C_d for Y-direction
rho_LE_x         = BVP.hyp.LE.rho(1); % ρ for X-direction
rho_LE_y         = BVP.hyp.LE.rho(2); % ρ for Y-direction
mean_u_active_LE = BVP.lsNIPCE.LE.mean_u_active_pc;
cov_u_active_LE  = BVP.lsNIPCE.LE.cov_u_active_pc;

% St. Venant (ST) model parameters
C_e_ST = BVP.hyp.ST.C_e;
C_d_ST_x = BVP.hyp.ST.C_d{1}; % C_d for X-direction
C_d_ST_y = BVP.hyp.ST.C_d{2}; % C_d for Y-direction
rho_ST_x = BVP.hyp.ST.rho(1); % ρ for X-direction
rho_ST_y = BVP.hyp.ST.rho(2); % ρ for Y-direction
mean_u_active_ST = BVP.lsNIPCE.ST.mean_u_active_pc;
cov_u_active_ST = BVP.lsNIPCE.ST.cov_u_active_pc;

%% Calculate Log-Likelihood for LE
funcTheta_LE = 0;

for iDof = 1:2
  % Assign direction-specific parameters
  if iDof == 1
    C_d_LE = C_d_LE_x;
    rho_LE = rho_LE_x;
  else
    C_d_LE = C_d_LE_y;
    rho_LE = rho_LE_y;
  end
  
  % Projection matrix for current DOF
  P = P_active{iDof};
  % Covariance matrix for LE
  Sigma_LE = rho_LE ^ 2 * P * cov_u_active_LE * P' + C_e_LE + C_d_LE;
  % Cholesky decomposition
  L_LE = chol(Sigma_LE, 'lower');
  inv_Sigma_LE = chol_solve(L_LE, eye(nSen));
  
  % Log-likelihood components
  alpha1_LE = nRep * nSen * log(2 * pi);
  alpha2_LE = 2 * nRep * sum(log(diag(L_LE)));
  alpha3_LE = 0;
  
  % Sum over repetitions
  for iRep = 1:nRep
    residual = Y_exp(iDof:2:end, iRep) - rho_LE * P * mean_u_active_LE;
    alpha3_LE = alpha3_LE + residual' * inv_Sigma_LE * residual;
  end
  
  funcTheta_LE = funcTheta_LE + 0.5 * (alpha1_LE + alpha2_LE + alpha3_LE);
end

%% Calculate Log-Likelihood for ST
funcTheta_ST = 0;

for iDof = 1:2
  % Assign direction-specific parameters
  if iDof == 1
    C_d_ST = C_d_ST_x;
    rho_ST = rho_ST_x;
  else
    C_d_ST = C_d_ST_y;
    rho_ST = rho_ST_y;
  end
  
  % Projection matrix for current DOF
  P = P_active{iDof};
  % Covariance matrix for ST
  Sigma_ST = rho_ST ^ 2 * P * cov_u_active_ST * P' + C_e_ST + C_d_ST;
  % Cholesky decomposition
  L_ST = chol(Sigma_ST, 'lower');
  inv_Sigma_ST = chol_solve(L_ST, eye(nSen));
  
  % Log-likelihood components
  alpha1_ST = nRep * nSen * log(2 * pi);
  alpha2_ST = 2 * nRep * sum(log(diag(L_ST)));
  alpha3_ST = 0;
  
  % Sum over repetitions
  for iRep = 1:nRep
    residual = Y_exp(iDof:2:end, iRep) - rho_ST * P * mean_u_active_ST;
    alpha3_ST = alpha3_ST + residual' * inv_Sigma_ST * residual;
  end
  
  funcTheta_ST = funcTheta_ST + 0.5 * (alpha1_ST + alpha2_ST + alpha3_ST);
end

% Convert to negative log-likelihood
funcTheta_LE = -funcTheta_LE;
funcTheta_ST = -funcTheta_ST;

%% Calculate Bayes Factor
% Logarithmic difference
differenceLog = (funcTheta_ST - funcTheta_LE) / (nRep * nSen);
% Bayes Factor
BF = exp(differenceLog);

%% Assign results back to BVP
BVP.modelSelection.BF = BF;
BVP.modelSelection.funcTheta_ST = funcTheta_ST;
BVP.modelSelection.funcTheta_LE = funcTheta_LE;

%% Display Results
fprintf('Log-Likelihood (ST): %f\n', funcTheta_ST);
fprintf('Log-Likelihood (LE): %f\n', funcTheta_LE);
fprintf('Bayes Factor (BF): %f\n', BF);

end
