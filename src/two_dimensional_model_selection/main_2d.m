
% MAIN SCRIPT FOR INFINITE PLATE BEAM ANALYSIS
% -------------------------------------------------------------------------
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
%
% Description:
% This script performs structural analysis on an infinite plate beam using:
%   - Deterministic FEM (Linear Elasticity and St. Venant-Kirchhoff models)
%   - Monte Carlo Simulations for uncertainty quantification
%   - Polynomial Chaos Expansion (PCE) via Spectral FEM
%   - Bayesian inference through statFEM with hyperparameter estimation
%   - Model comparison using RMSE and Bayes Factor calculations
%
% Features:
%   - Linear Elastic and St. Venant-Kirchhoff material models
%   - Monte Carlo simulations for stochastic analysis
%   - Polynomial Chaos expansions for uncertainty quantification
%   - Hyperparameter estimation using Bayesian methods
%   - Visualization of displacement and stress fields
%
% Dependencies:
%   • Gmsh (for mesh generation)
%   • MATLAB toolboxes: Statistics and Machine Learning, Optimization
%   • User-defined libraries under '../../lib/'
%
% Usage:
%   Run this script directly to perform the full analysis pipeline.
%
% -------------------------------------------------------------------------
clearvars; close all; clc;

% Add library path for required functions
addpath('../../lib/')

% Initialize Boundary Value Problem (BVP) structure
BVP = [];

%% ---------------------- PREPROCESSING -------------------------------
% FEM Preprocessing: Generate mesh, define material properties, and setup BCs
BVP = infinitePlate_FEM_preprocess(BVP);

%% ---------- FEM PROCESSING & POSTPROCESSING (DETERMINISTIC) ----------
% Linear Elastic Model
BVP = infinitePlate_FEM_processLE(BVP);
BVP = infinitePlate_FEM_postprocessLE(BVP);
BVP = infinitePlate_balancedStatus_LE(BVP); % Check equilibrium

% St. Venant-Kirchhoff Model
BVP = infinitePlate_FEM_processST(BVP);
BVP = infinitePlate_FEM_postprocessST(BVP);
BVP = infinitePlate_balancedStatus_ST(BVP); % Check equilibrium

%% ---------------- MONTE CARLO SIMULATIONS (OPTIONAL) ----------------
% Linear Elastic Model - Monte Carlo Simulation
BVP = infinitePlate_FEM_processLE_monteCarlo(BVP);

% St. Venant-Kirchhoff Model - Monte Carlo Simulation
BVP = infinitePlate_FEM_processST_monteCarlo(BVP);

%% ---- SPECTRAL FEM PROCESSING & POSTPROCESSING (POLYNOMIAL CHAOS) ----
% Linear Elastic Model (PC Expansion)
BVP = infinitePlate_FEM_processLE_NIPC(BVP);

% St. Venant-Kirchhoff Model (PC Expansion)
BVP = infinitePlate_FEM_processST_NIPC(BVP);

%% ------------------- OBSERVATION DATA GENERATION ---------------------
% Synthetic experimental data based on the St. Venant-Kirchhoff model
BVP = infinitePlate_obs_generate(BVP);

%% ------------------- PROJECTION MATRIX GENERATION --------------------
% Calculate the observation projection matrix (P)
BVP = infinitePlate_P_matrix(BVP);

%% --------------------- HYPERPARAMETER ESTIMATION ---------------------
% Estimate hyperparameters for statFEM (Bayesian framework)
BVP = infinitePlate_hyperParameter_ST(BVP); % For ST model
BVP = infinitePlate_hyperParameter_LE(BVP); % For LE model

%% ------ STATFEM PROCESSING & POSTPROCESSING (BAYESIAN INFERENCE) ------
% Linear Elastic Model with statFEM
BVP = infinitePlate_statFEM_processLE(BVP, 50);
BVP = infinitePlate_statFEM_postprocessLE(BVP);

% St. Venant-Kirchhoff Model with statFEM
BVP = infinitePlate_statFEM_processST(BVP, 50);
BVP = infinitePlate_statFEM_postprocessST(BVP);

% Optional: Check balanced status for statFEM solution
BVP = infinitePlate_statFEMbalancedStatus_ST(BVP);

%% ------------------------ MODEL EVALUATION ---------------------------
% Compute Root Mean Square Error (RMSE)
BVP = infinitePlate_RMSE(BVP);

% Calculate Bayes Factor for model comparison
BVP = infinitePlate_BayesFactor(BVP);

%% --------------------------- VISUALIZATION ---------------------------
% Generate surface plots for displacement and stress fields
BVP = infinitePlate_FEM_surfplot(BVP, 1);

% ----------------------------- END ------------------------------------
disp('Infinite Plate Beam Analysis Completed Successfully!');

