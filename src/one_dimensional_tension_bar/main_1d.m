
% main_1D.m
% -------------------------------------------------------------------------
% Main script for the simulation of a 1D tension bar under a tip load using
% both deterministic FEM, Monte Carlo simulations, Polynomial Chaos (PC),
% and Statistical FEM (StatFEM) approaches.
%
% The script performs the following tasks:
%   - Preprocessing: Mesh generation and problem setup
%   - FEM Processing: Deterministic linear elasticity analysis
%   - Monte Carlo Simulations: Propagation of uncertainty via random samples
%   - Polynomial Chaos Expansion: Non-intrusive PCE for uncertainty quantification
%   - Experimental Data Generation: Synthetic observations (linear/nonlinear)
%   - Projection Matrix Calculation: Constructs projection matrix for sensor data
%   - Hyperparameter Estimation: Estimates noise/model uncertainty parameters
%   - StatFEM Processing: Assimilates experimental data with FEM predictions
%   - Results Visualization: Plots displacements and uncertainties
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% -------------------------------------------------------------------------

clearvars; close all; clc;

%% Add library path
% Add custom libraries required for functions used in the simulation
addpath('../../lib/')

%% Initialize Boundary Value Problem (BVP) structure
BVP = [];

%% Preprocessing
% Generate the mesh, define boundary conditions, and set material parameters.
BVP = tensionBar_1D_preprocess(BVP);

%% FEM Processing (Deterministic)
% Solve the boundary value problem using deterministic FEM.
BVP = tensionBar_1D_FEM_processLE(BVP);

%% FEM Monte Carlo Simulations (Linear Elasticity)
% Perform Monte Carlo analysis with randomly sampled Young's modulus.
BVP = tensionBar_1D_FEM_processLE_MC(BVP);

%% FEM PC Simulation (Linear Elasticity)
% Perform Polynomial Chaos Expansion (PCE) for uncertainty quantification.
BVP = tensionBar_1D_FEM_processLE_PC(BVP);

%% StatFEM Processing
% Generate synthetic experimental data (choose 'linear' or 'nonlinear').
% Uncomment the line below for linear data generation:
% BVP = tensionBar_1D_obs_generate(BVP, 'linear', 7);

% Generate nonlinear synthetic experimental data with specified case:
BVP = tensionBar_1D_obs_generate(BVP, 'nonlinear', 7);

%% Calculate Projection Matrix P
% Construct the projection matrix mapping FEM predictions to sensor locations.
BVP = tensionBar_1D_P_matrix(BVP);

%% Estimate hyperparameters for StatFEM
% Estimate noise, correlation length, and scaling factors from experimental data.
BVP = tensionBar_1D_hyperParameter_LE(BVP);

%% Perform StatFEM Analysis
% Assimilate experimental data with FEM predictions using identified hyperparameters.
BVP = tensionBar_1D_statFEM_processLE(BVP);

%% Plot Results
% Visualize the computed displacement fields, uncertainty bounds, and experimental data.
BVP = tensionBar_1D_plot(BVP);

%% Completion Message
disp('1D Tension Bar Simulation Completed!');

