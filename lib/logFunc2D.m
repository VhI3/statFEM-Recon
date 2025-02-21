function [f] = logFunc2D(params, Cu, P, C_e, y_i, U_mean, x, nSen, ii)
% LOGFUNC2D Computes the log-posterior objective function for 2D analysis.
%   This function evaluates the log-posterior objective function for a
%   2D problem, using input parameters, covariance matrices, and sensor data.
%
% Inputs:
%   params - Hyperparameter vector to evaluate.
%   Cu     - Covariance matrix if prior displacement.
%   P      - Projection matrix.
%   C_e    - Error covariance matrix.
%   y_i    - Sensor observations.
%   U_mean - Mean displacement field.
%   x      - Sensor locations.
%   nSen   - Number of sensors.
%   ii     - Index or additional input for specific functionality.
%
% Outputs:
%   f      - Computed log-posterior objective function value.
%
% Example Usage:
%   f = logFunc2D(params, Cu, P, C_e, y_i, U_mean, x, nSen, ii);
%
% Note: Gradient computation can be added if needed in future versions.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

% Calculate the objective function (log-posterior)
f = logposterior2D(params, Cu, P, C_e, y_i, U_mean, x, nSen, ii);

end
