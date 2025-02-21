function BVP = infinitePlate_P_matrix(BVP)
% INFINITEPLATE_P_MATRIX Constructs the observation matrix for FEM analysis.
%   This function creates a mapping matrix that associates sensor locations
%   with the corresponding FEM nodes and degrees of freedom (DOFs).
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       Sensor coordinates, node coordinates, DOFs, and active DOFs.
%
% Outputs:
%   BVP - Updated BVP structure with computed observation matrices.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
senCoor = BVP.obs.senCoor; % Sensor coordinates
nSen = BVP.obs.nSen; % Number of sensors
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
GDOFs = BVP.msh.GDOFs; % Global degrees of freedom
DOFs = BVP.msh.DOFs; % Degrees of freedom per node
activeDOFsX = BVP.BC.activeDOFsX; % Active DOFs in X-direction
activeDOFsY = BVP.BC.activeDOFsY; % Active DOFs in Y-direction
activeDOFs = BVP.BC.activeDOFs; % Active DOFs

%% Initialize Observation Matrix and Sensor Nodes
P_matrix = zeros(nSen * DOFs, GDOFs); % Observation matrix
senNode = zeros(nSen, 1); % Sensor-node mapping
cc = 1; % Counter for observation rows

%% Map Sensors to Nodes
for i = 1:nSen
  % Find nodes matching sensor coordinates (within tolerance)
  rows = find(abs(nodeCoordinates(:, 1) - senCoor(i, 1)) < eps);
  nodes = abs(nodeCoordinates(rows, 2) - senCoor(i, 2)) < eps;
  senNode(i) = rows(nodes);
  
  % Assign values in the observation matrix
  P_matrix(cc, DOFs * senNode(i) - 1) = 1; % X DOF
  P_matrix(cc + 1, DOFs * senNode(i)) = 1; % Y DOF
  cc = cc + 2; % Move to the next sensor
end

%% Create Active DOF Observation Matrices
Px = P_matrix(1:2:end, activeDOFsX); % X-component matrix
Py = P_matrix(2:2:end, activeDOFsY); % Y-component matrix
P_active{1} = P_matrix(1:2:end, activeDOFs); % Combined active X observations
P_active{2} = P_matrix(2:2:end, activeDOFs); % Combined active Y observations
P_matrix_active = P_matrix(:, activeDOFs); % Combined active observation matrix

%% Assign Results Back to BVP
BVP.obs.senNode = senNode; % Sensor-node mapping
BVP.obs.P_matrix = P_matrix; % Full observation matrix
BVP.obs.Px = Px; % X-component matrix
BVP.obs.Py = Py; % Y-component matrix
BVP.obs.P_active = P_active; % Active DOF observation matrices
BVP.obs.P_matrix_active = P_matrix_active; % Active DOF observation matrix

end
