function BVP = infinitePlate_FEM_processST_NIPC(BVP)
% INFINITEPLATE_FEM_PROCESSST_NIPC Processes FEM results for nonlinear analysis using NIPCE.
%   This function calculates the displacement response using the
%   non-intrusive polynomial chaos expansion (NIPCE) approach for
%   nonlinear steady-state analysis.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       lsNIPCE parameters, mesh details, boundary conditions, and material properties.
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
Xi_PC = BVP.lsNIPCE.Xi_PC; % Input random variables
E_PC = BVP.lsNIPCE.E_PC; % Young's modulus realizations

% Mesh and boundary conditions
GDOFs = BVP.msh.GDOFs; % Global degrees of freedom
activeDOFsX = BVP.BC.activeDOFsX; % Active DOFs in X direction
activeDOFsY = BVP.BC.activeDOFsY; % Active DOFs in Y direction
activeDOFs = BVP.BC.activeDOFs; % Active DOFs
bottomNodesDofx = BVP.msh.bottomNodesDofx; % Bottom node DOFs in X direction
elementNodes = BVP.msh.elementNodes; % Element connectivity

%% Initialize Storage for Displacement
bigPsi = zeros(N_s, P_PCE); % Basis matrix
u_realization = zeros(N_s, GDOFs); % Displacement realizations
cov_u_pc = zeros(GDOFs); % Covariance matrix for displacement

%% Process Using Least-Squares NIPCE
pb = CmdLineProgressBar('Processing NIPCE for St. Venant-Kirchhoff...');
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
  BVP_tmp = infinitePlate_FEM_processST(BVP_tmp);
  u_realization(i, :) = BVP_tmp.proc.ST.u;
end

%% Compute Statistical Quantities
% Mean displacement coefficients
u_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * u_realization);
u_NIPC = u_NIPC';
mean_u_pc = u_NIPC(:, 1);
mean_ux_pc = mean_u_pc(1:2:end);
mean_uy_pc = mean_u_pc(2:2:end);

% Covariance of displacement
for j = 2:P_PCE
  cov_u_pc = cov_u_pc + (1 / sqrt(factorial(j))) ^ 2 * PsiSqNorm(j) * u_NIPC(:, j) * u_NIPC(:, j)';
end

cov_u_pc = 0.5 * (cov_u_pc + cov_u_pc'); % Ensure symmetry
ci_cov_u_pc = sqrt(diag(cov_u_pc)) * 1.96;

% Surface representation
mean_ux_pc_Surf = makeSurf(elementNodes, mean_ux_pc);
std_ux_pc = sqrt(diag(cov_u_pc(1:2:end, 1:2:end)));
std_ux_pc_Surf = makeSurf(elementNodes, std_ux_pc);

%% Assign Results Back to BVP
BVP.lsNIPCE.ST.mean_u_pc = mean_u_pc;
BVP.lsNIPCE.ST.mean_ux_pc = mean_ux_pc;
BVP.lsNIPCE.ST.mean_uy_pc = mean_uy_pc;
BVP.lsNIPCE.ST.mean_ux_active_pc = mean_u_pc(activeDOFsX);
BVP.lsNIPCE.ST.mean_uy_active_pc = mean_u_pc(activeDOFsY);
BVP.lsNIPCE.ST.mean_u_active_pc = mean_u_pc(activeDOFs);
BVP.lsNIPCE.ST.mean_ux_bottomNode_pc = mean_u_pc(bottomNodesDofx);
BVP.lsNIPCE.ST.mean_ux_pc_Surf = mean_ux_pc_Surf;

BVP.lsNIPCE.ST.cov_u_pc = cov_u_pc;
BVP.lsNIPCE.ST.ci_cov_u_pc = ci_cov_u_pc;
BVP.lsNIPCE.ST.cov_ux_pc = cov_u_pc(1:2:end, 1:2:end);
BVP.lsNIPCE.ST.std_ux_pc = std_ux_pc;
BVP.lsNIPCE.ST.ci_cov_ux_pc = std_ux_pc * 1.96;
BVP.lsNIPCE.ST.cov_uy_pc = cov_u_pc(2:2:end, 2:2:end);
BVP.lsNIPCE.ST.std_uy_pc = sqrt(diag(cov_u_pc(2:2:end, 2:2:end)));
BVP.lsNIPCE.ST.ci_cov_uy_pc = sqrt(diag(cov_u_pc(2:2:end, 2:2:end))) * 1.96;
BVP.lsNIPCE.ST.cov_ux_active_pc = cov_u_pc(activeDOFsX, activeDOFsX);
BVP.lsNIPCE.ST.cov_uy_active_pc = cov_u_pc(activeDOFsY, activeDOFsY);
BVP.lsNIPCE.ST.cov_u_active_pc = cov_u_pc(activeDOFs, activeDOFs);
BVP.lsNIPCE.ST.ci_cov_ux_bottomNode_pc = ci_cov_u_pc(bottomNodesDofx);
BVP.lsNIPCE.ST.std_ux_pc_Surf = std_ux_pc_Surf;

end
