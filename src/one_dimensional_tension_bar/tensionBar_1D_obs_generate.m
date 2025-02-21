function BVP = tensionBar_1D_obs_generate(BVP, obsCase, calCase)
% TENSIONBAR_1D_OBS_GENERATE
% Generates synthetic experimental data for a 1D tension bar under tip load.
%
% This function simulates noisy experimental observations using a log-normal
% model with correlated noise and different sensor placements.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing problem parameters.
%
% Outputs:
%   BVP - Updated BVP structure with generated observations.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)



%% Step 1: Define Experiment Parameters
nrep_sample = 1000; % Number of repetitions (noise mismatched realizations)
n_e_sample = 100; % Number of evaluation points

L = BVP.geometry.L;
A = BVP.geometry.A;
mu_E = BVP.material.mu_E;
f_bar = BVP.loading.f_bar;

xx = linspace(0, L, n_e_sample + 1)'; % Discretized positions

%% Step 2: Define True Displacement Model
if strcmp(obsCase, 'linear')
    S = 0;
    u_sample = (f_bar / (A * mu_E)) * xx;
elseif strcmp(obsCase, 'nonlinear')
    S = 0.015; % Strain scaling factor
    u_sample = (f_bar / (A * S * mu_E)) * (1 - exp(-S * xx));
end

%% Step 3: Noise Model & Hyperparameters
noiseVector = [0.004; 0.04; 0.4; 4]; % Noise levels
epsi = noiseVector(1); % Select noise level
rho_sample = 1.2; % Scaling factor
sig_d_sample = 0.9; % Variance of model-reality mismatch
l_d_sample = 4; % Correlation length of model-reality mismatch

%% Step 4: Generate Experimental Noise
rng(4, 'twister'); % Set random seed for reproducibility
C_e_sample = epsi * eye(length(xx)); % Observation noise covariance
e_sample = mvnrnd(zeros(length(xx), 1), C_e_sample, nrep_sample)';

rng(5, 'twister');
C_d_sample = sqexp(xx, xx, log(sig_d_sample), log(l_d_sample)); % Correlated noise
d_sample = mvnrnd(zeros(length(xx), 1), C_d_sample, nrep_sample)';

%% Step 5: Generate Experimental Observations
Y_experiment = zeros(length(xx), nrep_sample);

for i = 1:nrep_sample
    Y_experiment(:, i) = rho_sample * u_sample + d_sample(:, i) + e_sample(:, i);
end

%% Step 6: Select Sensor Positions Based on Experiment Case

% Number of repetitions and sensors
switch calCase
    case 1, nrep = 1; nsen = 11;
    case 2, nrep = 10; nsen = 11;
    case 3, nrep = 100; nsen = 11;
  case 4, nrep = 1000; nsen = 11;
  case 5, nrep = 1; nsen = 33;
  case 6, nrep = 10; nsen = 33;
  case 7, nrep = 100; nsen = 33;
  case 8, nrep = 1000; nsen = 33;
  case 9, nrep = 1; nsen = 50;
  case 10, nrep = 10; nsen = 50;
  case 11, nrep = 100; nsen = 50;
  case 12, nrep = 1000; nsen = 50;
  otherwise , error('Invalid calCase selected.');
end

% Define sensor placement strategy
if nsen == 4
  senInd = [20, 40, 60, 80];
elseif nsen == 11
  senInd = [5, 20, 25, 35, 40, 50, 60, 75, 80, 90, 101];
elseif nsen == 33
  senInd = round(linspace(2, 101, nsen));
else
  senInd = 1:2:100;
end

% Extract selected sensor positions and corresponding observations
senCoor = xx(senInd);
y_obs = Y_experiment(senInd, 1:nrep);

C_e = epsi * eye(length(senCoor), length(senCoor));

%% **Step 7: Assign Data to BVP Structure**
BVP.obs.nrep = nrep; % Number of repetitions used
BVP.obs.nsen = nsen; % Number of sensors used
BVP.obs.senInd = senInd; % Sensor indices
BVP.obs.senCoor = senCoor; % Sensor coordinates
BVP.obs.Y_experiment = Y_experiment; % Full dataset
BVP.obs.y_obs = y_obs; % Selected observations
BVP.obs.epsi = epsi; % Observation noise level
BVP.obs.rho_sample = rho_sample; % Scaling factor
BVP.obs.sig_d_sample = sig_d_sample; % Correlated noise std dev
BVP.obs.l_d_sample = l_d_sample; % Correlation length scale
BVP.obs.C_e_sample = C_e_sample; % Observation noise covariance
BVP.obs.C_d_sample = C_d_sample; % Correlated noise covariance
BVP.obs.C_e = C_e; % Covariance of Observation noise
BVP.obs.S = S; % Strain scaling factor

disp('5. Synthetic Observations Generated.');
end
