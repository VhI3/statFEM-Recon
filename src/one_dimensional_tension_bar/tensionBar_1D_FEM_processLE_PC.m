function BVP = tensionBar_1D_FEM_processLE_PC(BVP)
% tensionBar_1D_FEM_processLE_PC
% Performs Polynomial Chaos Expansion (PCE) for a 1D tension bar under a tip load.
% This function uses a non-intrusive PCE approach to quantify uncertainty
% in displacement due to variability in Young’s modulus.
%
% Key Steps:
%   1. Generate Hermite polynomial chaos basis functions.
%   2. Perform FEM simulations for different realizations of Young's modulus.
%   3. Estimate PCE coefficients using least squares regression.
%   4. Compute the mean and covariance of the displacement field.
%   5. Evaluate PCE for generating displacement samples.
%   6. Estimate PDF and CDF for the displacement at the tip.
%
% Inputs:
%   BVP - Structure containing the boundary value problem setup and parameters.
%
% Outputs:
%   BVP - Updated structure with PCE results and uncertainty quantification.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign Parameters
nPC = BVP.UQ.nPC;
nPC_s = BVP.UQ.nPC_s;
P_PCE = BVP.UQ.P_PCE;
Xi_PC = BVP.UQ.Xi_PC;
E_PC = BVP.UQ.E_PC;
numberNodes = BVP.mesh.numberNodes;
f_bar = BVP.loading.f_bar;
A = BVP.geometry.A;
nodeCoordinates = BVP.mesh.nodeCoordinates;
nElm = BVP.mesh.nElm;
elementNodes = BVP.fem.elementNodes; % Element connectivity
activeDOFs = BVP.fem.activeDOFs; % Active degrees of freedom
GDof = BVP.fem.GDof; % Global degrees of freedom
%%
% Generate Hermite Polynomial Chaos Basis
[~, Psi_s, ~, PsiSqNorm, P_PCE] = Hermite_PC(1, P_PCE);
bigPsi = zeros(nPC, P_PCE);

for i = 1:nPC
    
    for j = 1:P_PCE
        bigPsi(i, j) = double(subs(Psi_s{j, 1}, Xi_PC(i))); % Evaluate Hermite polynomials
    end
    
end

% FEM Analysis for Each Sample
u_k = zeros(nPC, numberNodes);

for i = 1:nPC
    u_k(i, :) = FEM_Bar_deter_Tipload(E_PC(i), f_bar, A, nElm, elementNodes, nodeCoordinates, activeDOFs, GDof);
end

% Solve for PCE Coefficients using Least Squares Regression
u_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * u_k);
u_NIPC = u_NIPC';

% Mean Displacement from PC Expansion
mu_u_pc = u_NIPC(:, 1);
mu_u_active_pc = mu_u_pc(activeDOFs);

% Compute Covariance Matrix
C_u_NIPC = zeros(BVP.mesh.numberNodes);

for j = 2:P_PCE
    C_u_NIPC = C_u_NIPC + PsiSqNorm(j) * u_NIPC(:, j) * u_NIPC(:, j)';
end

C_u_pc = 0.5 * (C_u_NIPC + C_u_NIPC');
C_u_active_pc = C_u_pc(activeDOFs, activeDOFs);
ci_u_pc = sqrt(diag(C_u_pc)) * 1.96;

%% Representation of the Response Displacement
% Initialize displacement storage
u_samples = zeros(numberNodes, nPC_s);

% Generate standard normal samples (Gaussian RV)
xi_samples = randn(nPC_s, 1);

% Define symbolic variable for PCE evaluation
xi_v = sym('xi_1'); % Symbolic variable
xi_n = xi_samples; % Numeric samples

% Evaluate PCE Expansion for Displacement
for j = 0:(P_PCE - 1)
    fprintf('Evaluating PC expansion %g/%g\n', j, P_PCE - 1);
    
    % Extract PCE polynomial basis function
    psi_j = double(subs(Psi_s{j + 1}, xi_v, xi_n));
    
    % Extract corresponding PCE coefficient
    d_j = u_NIPC(:, j + 1);
    
    % Compute displacement realization
    for i = 1:numberNodes
        u_samples(i, :) = u_samples(i, :) + d_j(i) * psi_j';
    end
    
end

%% **Extract & Analyze Displacement at the End Node**
u_PC = u_samples(end, :); % Displacement at the last node

% Compute Probability Density Function (PDF)
[pdf_u_tip_pc, xpc_pdf] = ksdensity(u_PC); % pdf

% Compute Cumulative Distribution Function (CDF)
[cdf_u_tip_pc, xpc_cdf] = ecdf(u_PC); % cdf
%% Assign Results to BVP Structure
BVP.UQ.P_PCE = P_PCE;
BVP.UQ.nPC = nPC;
BVP.UQ.Xi_PC = Xi_PC;
BVP.UQ.E_PC = E_PC;
BVP.UQ.bigPsi = bigPsi;
BVP.UQ.u_NIPC = u_NIPC;
BVP.UQ.mu_u_pc = mu_u_pc;
BVP.UQ.C_u_pc = C_u_pc;
BVP.UQ.ci_u_pc = ci_u_pc;
BVP.UQ.LE.pdf_u_tip_pc = pdf_u_tip_pc;
BVP.UQ.LE.cdf_u_tip_pc = cdf_u_tip_pc;
BVP.UQ.LE.xpc_pdf = xpc_pdf;
BVP.UQ.LE.xpc_cdf = xpc_cdf;
BVP.UQ.LE.mu_u_active_pc = mu_u_active_pc;
BVP.UQ.LE.C_u_active_pc = C_u_active_pc;
disp('4. Polynomial Chaos Simulation for Linear Elasticity Completed.');
end
