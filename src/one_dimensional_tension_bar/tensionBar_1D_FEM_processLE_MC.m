function BVP = tensionBar_1D_FEM_processLE_MC(BVP)
% Monte Carlo Simulation for 1D Tension Bar with Linear Elasticity
%
% This function performs Monte Carlo simulations for a 1D tension bar under
% a tip load by varying Young's modulus (`E`) based on stochastic realizations.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP (Extract Input Parameters)
nMC = BVP.UQ.nMC; % Number of Monte Carlo Samples
E_MC = BVP.UQ.E_MC; % Monte Carlo realizations of Young’s modulus
f_bar = BVP.loading.f_bar; % Tip load applied at the free end
A = BVP.geometry.A; % Cross-sectional area
nElm = BVP.mesh.nElm; % Number of elements
numberNodes = BVP.mesh.numberNodes; % Total number of nodes
nodeCoordinates = BVP.mesh.nodeCoordinates; % Node positions along the bar
elementNodes = BVP.fem.elementNodes; % Element connectivity
activeDOFs = BVP.fem.activeDOFs; % Active degrees of freedom
GDof = BVP.fem.GDof; % Global degrees of freedom

%% Monte Carlo Simulation Loop
pb = CmdLineProgressBar('Monte Carlo Simulation for Linear Elasticity...');
u_mc = zeros(nMC, numberNodes); % Store displacements for all MC samples

for i = 1:nMC
    pb.print(i, nMC); % Display progress
    
    % Call FEM solver with a different Young's modulus for each realization
    u_mc(i, :) = FEM_Bar_deter_Tipload(E_MC(i), f_bar, A, nElm, elementNodes, nodeCoordinates, activeDOFs, GDof);
end

%% Compute Statistical Measures
mean_u_mc = mean(u_mc); % Mean displacement
std_u_mc = std(u_mc); % Standard deviation of displacement

[pdf_u_tip_mc, xmc_pdf] = ksdensity(u_mc(:, end)); % pdf
[cdf_u_tip_mc, xmc_cdf] = ecdf(u_mc(:, end)); % cdf

%% Assign Results Back to BVP
BVP.UQ.LE.u_mc = u_mc;
BVP.UQ.LE.mean_u_mc = mean_u_mc;
BVP.UQ.LE.std_u_mc = std_u_mc;
BVP.UQ.LE.pdf_u_tip_mc = pdf_u_tip_mc;
BVP.UQ.LE.cdf_u_tip_mc = cdf_u_tip_mc;
BVP.UQ.LE.xmc_pdf = xmc_pdf;
BVP.UQ.LE.xmc_cdf = xmc_cdf;

disp('3. Monte Carlo Simulation for Linear Elasticity Completed.');
end
