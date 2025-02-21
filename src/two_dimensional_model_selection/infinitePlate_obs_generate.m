function BVP = infinitePlate_obs_generate(BVP)
% INFINITEPLATE_OBS_GENERATE Generates synthetic observations for FEM analysis.
%   This function generates synthetic sensor observations based on FEM
%   results, incorporating sampling error and noise.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       FEM results (displacements), observation parameters, and sensor details.
%
% Outputs:
%   BVP - Updated BVP structure with synthetic observation data.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
ux = BVP.proc.ST.ux; % X-direction displacement
uy = BVP.proc.ST.uy; % Y-direction displacement
nSam = BVP.obs.nSam; % Number of samples
d_sample = BVP.obs.d_sample; % Sampling error
e_sample = BVP.obs.e_sample; % Observation noise
rho_sample = BVP.obs.rho_sample; % Scaling factor
nSen = BVP.obs.nSen; % Number of sensors
senIndex = BVP.obs.sensorIndex; % Sensor indices
bottomSen_DofX = BVP.obs.bottomSen_DofX; % Sensor indices for bottom nodes (X)

%% Initialize Synthetic Observations
Y_exp_x = zeros(nSen, nSam); % Observations in X-direction
Y_exp_y = zeros(nSen, nSam); % Observations in Y-direction

%% Generate Observations
for i = 1:nSam
    % Synthetic observations with sampling error and noise
    Y_exp_x(:, i) = rho_sample * ux(senIndex) + d_sample(:, i) + e_sample(:, i);
    Y_exp_y(:, i) = rho_sample * uy(senIndex) + d_sample(:, i) + e_sample(:, i);
end

% Combine X and Y observations into a single array
Y_exp = zeros(2 * nSen, nSam);
Y_exp(1:2:end, :) = Y_exp_x; % Fill X observations
Y_exp(2:2:end, :) = Y_exp_y; % Fill Y observations

% Extract observations for bottom nodes (X-direction)
Y_exp_bottom_DofX = Y_exp(bottomSen_DofX, :);

% Compute mean of bottom node observations
mean_Y_exp_bottom_DofX = mean(Y_exp_bottom_DofX, 2);

%% Assign Results Back to BVP
BVP.obs.Y_exp_x = Y_exp_x; % Observations in X-direction
BVP.obs.Y_exp_y = Y_exp_y; % Observations in Y-direction
BVP.obs.Y_exp = Y_exp; % Combined observations
BVP.obs.Y_exp_bottom_DofX = Y_exp_bottom_DofX; % Bottom node observations
BVP.obs.mean_Y_exp_bottom_DofX = mean_Y_exp_bottom_DofX; % Mean of bottom node observations

end
