function BVP = infinitePlate_FEM_processLE_NIPC(BVP)
% INFINITEPLATE_FEM_PROCESSLE_NIPC Processes FEM results using non-intrusive PCE.
%   This function implements a non-intrusive polynomial chaos expansion (NIPCE)
%   approach to calculate the displacement and strain response of an infinite plate
%   under linear elasticity.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       lsNIPCE parameters, mesh details, boundary conditions, material properties.
%
% Outputs:
%   BVP - Updated BVP structure with computed response statistics.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
% NIPCE parameters
N_s = BVP.lsNIPCE.N_s; % Number of samples
P_PCE = BVP.lsNIPCE.P_PCE; % Number of PCE terms
Psi_s = BVP.lsNIPCE.Psi_s; % Polynomial chaos basis
PsiSqNorm = BVP.lsNIPCE.PsiSqNorm; % Norm of basis polynomials
Xi_PC = BVP.lsNIPCE.Xi_PC; % Input random variables (PCE samples)
E_PC = BVP.lsNIPCE.E_PC; % Young's modulus realizations

% Mesh and boundary conditions
GDOFs = BVP.msh.GDOFs; % Global degrees of freedom
activeDOFsX = BVP.BC.activeDOFsX; % Active DOFs in X direction
activeDOFsY = BVP.BC.activeDOFsY; % Active DOFs in Y direction
activeDOFs = BVP.BC.activeDOFs; % Active DOFs
elementNodes = BVP.msh.elementNodes; % Element connectivity

%% Initialize Storage for Displacement and Strain
bigPsi = zeros(N_s, P_PCE); % Basis matrix
u_realization = zeros(N_s, GDOFs); % Displacement realizations
strain_xx_realization = zeros(N_s, GDOFs / 2); % Strain in xx direction
strain_yy_realization = zeros(N_s, GDOFs / 2); % Strain in yy direction
cov_u_pc = zeros(GDOFs); % Covariance matrix for displacement

%% Process Using Least-Squares NIPCE
pb = CmdLineProgressBar('Processing PC for Linear Elasticity...');
BVP_tmp = BVP;

for i = 1:N_s
  pb.print(i, N_s);
  
  % Compute basis functions for the sample
  for j = 1:P_PCE
    Xi_PC_sym = sym(sprintf('%.17g', Xi_PC(i))); % Convert to a symbolic rational value
    bigPsi(i, j) = (1 / sqrt(factorial(j))) * double(subs(Psi_s{j, 1}, Xi_PC_sym));
    
    %bigPsi(i, j) = (1 / sqrt(factorial(j))) * double(subs(Psi_s{j, 1}, Xi_PC(i)));
  end
  
  % Update Young's modulus for the realization
  BVP_tmp.material.E = E_PC(i);
  
  % Solve FEM problem for the current realization
  BVP_tmp = infinitePlate_FEM_processLE(BVP_tmp);
  u_realization(i, :) = BVP_tmp.proc.LE.u;
  
  % Post-process for strain
  BVP_tmp = infinitePlate_FEM_postprocessLE(BVP_tmp);
  strain_xx_realization(i, :) = BVP_tmp.postproc.LE.strain_xx;
  strain_yy_realization(i, :) = BVP_tmp.postproc.LE.strain_yy;
end

%% Compute Statistical Quantities
% Mean displacement coefficients
u_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * u_realization);
u_NIPC = u_NIPC';
mean_u_pc = u_NIPC(:, 1);
mean_ux_pc = mean_u_pc(1:2:end);
mean_uy_pc = mean_u_pc(2:2:end);
mean_ux_active_pc = mean_u_pc(activeDOFsX);
mean_uy_active_pc = mean_u_pc(activeDOFsY);
mean_u_active_pc = mean_u_pc(activeDOFs);
% Covariance of displacement
for j = 2:P_PCE
  cov_u_pc = cov_u_pc + (1 / sqrt(factorial(j))) ^ 2 * PsiSqNorm(j) * u_NIPC(:, j) * u_NIPC(:, j)';
end

cov_u_pc = 0.5 * (cov_u_pc + cov_u_pc'); % Ensure symmetry
ci_cov_u_pc = sqrt(diag(cov_u_pc)) * 1.96;
cov_ux_active_pc = cov_u_pc(activeDOFsX, activeDOFsX);
cov_uy_active_pc = cov_u_pc(activeDOFsY, activeDOFsY);
cov_u_active_pc = cov_u_pc(activeDOFs, activeDOFs);
% Surface representation
mean_ux_pc_Surf = makeSurf(elementNodes, mean_ux_pc);
std_ux_pc = sqrt(diag(cov_u_pc(1:2:end, 1:2:end)));
std_ux_pc_Surf = makeSurf(elementNodes, std_ux_pc);

%% Strain Statistics
strain_xx_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * strain_xx_realization);
strain_NIPC_xx = strain_xx_NIPC';
mean_strain_xx = strain_NIPC_xx(:, 1);

strain_yy_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * strain_yy_realization);
strain_NIPC_yy = strain_yy_NIPC';
mean_strain_yy = strain_NIPC_yy(:, 1);

%% Assign Results Back to BVP
BVP.lsNIPCE.LE.mean_u_pc = mean_u_pc;
BVP.lsNIPCE.LE.mean_ux_pc = mean_ux_pc;
BVP.lsNIPCE.LE.mean_uy_pc = mean_uy_pc;
BVP.lsNIPCE.LE.mean_ux_pc_Surf = mean_ux_pc_Surf;
BVP.lsNIPCE.LE.mean_u_active_pc = mean_u_active_pc;
BVP.lsNIPCE.LE.mean_ux_active_pc = mean_ux_active_pc;
BVP.lsNIPCE.LE.mean_uy_active_pc = mean_uy_active_pc;
BVP.lsNIPCE.LE.std_ux_pc_Surf = std_ux_pc_Surf;
BVP.lsNIPCE.LE.strain_xx_NIPC = strain_xx_NIPC;
BVP.lsNIPCE.LE.mean_strain_xx = mean_strain_xx;
BVP.lsNIPCE.LE.strain_yy_NIPC = strain_yy_NIPC;
BVP.lsNIPCE.LE.mean_strain_yy = mean_strain_yy;
BVP.lsNIPCE.LE.cov_u_pc = cov_u_pc;
BVP.lsNIPCE.LE.ci_cov_u_pc = ci_cov_u_pc;
BVP.lsNIPCE.LE.cov_ux_active_pc = cov_ux_active_pc;
BVP.lsNIPCE.LE.cov_uy_active_pc = cov_uy_active_pc;
BVP.lsNIPCE.LE.cov_u_active_pc = cov_u_active_pc;
end
