function BVP = tensionBar_1D_preprocess(BVP)
% tensionBar_1D_preprocess
% Preprocessing function for the 1D tension bar problem under tip load.
%
% This function sets up the geometry, material properties, finite element
% mesh, and uncertainty quantification (UQ) parameters, including Monte Carlo
% simulations and Polynomial Chaos Expansion (PCE) settings.
%
% Key Features:
%   - Defines geometry, loading, and material parameters.
%   - Generates finite element mesh and degrees of freedom.
%   - Prepares Monte Carlo samples of Young’s modulus.
%   - Prepares samples for Polynomial Chaos Expansion.
%
% Inputs:
%   BVP - Boundary Value Problem structure (can be empty for initialization).
%
% Outputs:
%   BVP - Updated BVP structure containing:
%         * Geometry and loading conditions
%         * Material properties
%         * Mesh and FEM parameters
%         * Monte Carlo and PCE samples for UQ
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon


%% FEM Parameters
L = 100; % Length of the domain [mm]
f_bar = 800; % Magnitude of distributed force [N]
A = 20; % Cross-section area [mm^2]
nElm = 30; % Number of elements
nenod = 2; % Nodes per element
nodeCoordinates = linspace(0, L, nElm + 1)';
numberNodes = nElm + 1;
mu_E = 200; % Mean value of Young's modulus [MPa]
sig_E = 15; % Standard deviation of Young's modulus
exactSolution = (1 / (mu_E * A)) * (f_bar * nodeCoordinates);

DOFs = 1; % Degree of Freedom per node
GDof = DOFs * numberNodes; % Global degrees of freedom
fixedNode = 1;
activeDOFs = setdiff(1:GDof, fixedNode);
elementNodes = [(1:nElm)', (2:nElm + 1)'];

% Generate Monte Carlo Realizations of Young’s Modulus (E_MC)
% Number of Monte Carlo samples
nMC = 2000;

% Convert to Log-Normal Distribution Parameters
lambda = log(mu_E ^ 2 / sqrt(mu_E ^ 2 + sig_E ^ 2));
zeta = sqrt(log(1 + sig_E ^ 2 / mu_E ^ 2));

% Set random seed for reproducibility
rng(3, 'twister');

% Generate uniform random samples
uniformSample = rand(nMC, 1);

% Transform uniform samples to standard normal distribution
Xi_MC = norminv(uniformSample, 0, 1);

% Generate Young’s modulus realizations from log-normal distribution
E_MC = exp(lambda + zeta .* Xi_MC);

%% PCE Parameters
P_PCE = 8; % Order of Polynomial Chaos Expansion
nPC = 3 * P_PCE; % Number of samples for PCE
nPC_s = 2000; % Number of samples from PCE out

% Generate PCE Samples
rng(3, 'twister');
uniformSample = rand(nPC, 1);
Xi_PC = norminv(uniformSample, 0, 1); % Standard normal samples
E_PC = exp(lambda + zeta .* Xi_PC); % Log-normal realizations

%% Assign the BVP structure to the output variable
BVP.geometry.L = L;
BVP.loading.f_bar = f_bar;
BVP.geometry.A = A;
BVP.mesh.nElm = nElm;
BVP.mesh.nenod = nenod;
BVP.mesh.nodeCoordinates = nodeCoordinates;
BVP.mesh.numberNodes = numberNodes;
BVP.material.mu_E = mu_E;
BVP.material.sig_E = sig_E;
BVP.fem.activeDOFs = activeDOFs;
BVP.fem.GDof = GDof;
BVP.fem.DOFs = DOFs;
BVP.fem.elementNodes = elementNodes;
BVP.exactSolution = exactSolution;
BVP.UQ.nMC = nMC;
BVP.UQ.E_MC = E_MC;
BVP.UQ.Xi_MC = Xi_MC;
BVP.UQ.P_PCE = P_PCE;
BVP.UQ.nPC = nPC;
BVP.UQ.E_PC = E_PC;
BVP.UQ.Xi_PC = Xi_PC;
BVP.UQ.nPC_s = nPC_s;
%
disp('1. Preprocessing Completed')
end
