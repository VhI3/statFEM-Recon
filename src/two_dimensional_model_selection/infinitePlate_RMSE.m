function BVP = infinitePlate_RMSE(BVP)
% INFINITEPLATE_RMSE Calculates the Root Mean Square Error (RMSE) for LE and St. Venant (ST) models.
%
% This function computes the RMSE between the predicted displacements
% from statFEM with the St. Venant material model (ST) and linear elasticity (LE)
% models, and the observed data.
%
% Inputs:
%   BVP - Boundary value problem structure with model predictions and observations.
%
% Outputs:
%   BVP - Updated BVP structure containing RMSE results for ST and LE models.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
mean_UZx_y_ST = BVP.statFEM.ST.mean_UZx_y_red; % Predicted X-displacements (St. Venant)
mean_UZy_y_ST = BVP.statFEM.ST.mean_UZy_y_red; % Predicted Y-displacements (St. Venant)
mean_UZx_y_LE = BVP.statFEM.LE.mean_UZx_y_red; % Predicted X-displacements (Linear Elasticity)
mean_UZy_y_LE = BVP.statFEM.LE.mean_UZy_y_red; % Predicted Y-displacements (Linear Elasticity)
Y_exp_x = BVP.obs.Y_exp_x; % Observed X-displacements
Y_exp_y = BVP.obs.Y_exp_y; % Observed Y-displacements
nRep = BVP.statFEM.ST.nAssimilation; % Number of assimilation steps
nSen = BVP.obs.nSen; % Number of sensors

%% Calculate RMSE
summ_ST = 0;
summ_LE = 0;

for i = 1:nRep
    % Compute RMSE for St. Venant (ST) model
    rmse_ST_x = norm(mean_UZx_y_ST - Y_exp_x(:, i)) ^ 2; % X-component error
    rmse_ST_y = norm(mean_UZy_y_ST - Y_exp_y(:, i)) ^ 2; % Y-component error
    summ_ST = summ_ST + sqrt((rmse_ST_x + rmse_ST_y) / nSen);
    
    % Compute RMSE for Linear Elasticity (LE) model
    rmse_LE_x = norm(mean_UZx_y_LE - Y_exp_x(:, i)) ^ 2; % X-component error
    rmse_LE_y = norm(mean_UZy_y_LE - Y_exp_y(:, i)) ^ 2; % Y-component error
    summ_LE = summ_LE + sqrt((rmse_LE_x + rmse_LE_y) / nSen);
end

% Average RMSE across all assimilation steps
RMSE_ST = summ_ST / nRep;
RMSE_LE = summ_LE / nRep;

%% Display RMSE Results
fprintf('RMSE (St. Venant): %f\n', RMSE_ST);
fprintf('RMSE (Linear Elasticity): %f\n', RMSE_LE);

%% Assign back to BVP
BVP.RMSE.RMSE_ST = RMSE_ST;
BVP.RMSE.RMSE_LE = RMSE_LE;

end
