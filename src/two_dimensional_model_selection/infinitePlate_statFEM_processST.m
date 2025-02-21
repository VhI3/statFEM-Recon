function BVP = infinitePlate_statFEM_processST(BVP, n_0)
% INFINITEPLATE_STATFEM_PROCESSST Computes the statFEM posterior solution for St. Venant material.
%   This function computes the posterior mean and covariance for the
%   displacements in a St. Venant finite element model using statistical assimilation.
%
% Inputs:
%   BVP  - Boundary value problem structure.
%   n_0  - Number of reading from sensors.
%
% Outputs:
%   BVP  - Updated BVP structure containing posterior results.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
% Degrees of freedom
GDOFs = BVP.msh.GDOFs;
activeDOFsX = BVP.BC.activeDOFsX;
activeDOFsY = BVP.BC.activeDOFsY;
prescribedDOFsX = BVP.BC.prescribedDOFsX;
prescribedDOFsY = BVP.BC.prescribedDOFsY;

% Observation-related data
bottomNodesDofx = BVP.msh.bottomNodesDofx;
elementNodes = BVP.msh.elementNodes;
Px = BVP.obs.Px;
Py = BVP.obs.Py;

% Covariances and means from lsNIPCE
cov_ux_active = BVP.lsNIPCE.ST.cov_ux_active_pc;
cov_uy_active = BVP.lsNIPCE.ST.cov_uy_active_pc;
mean_ux_red = BVP.lsNIPCE.ST.mean_ux_active_pc;
mean_uy_red = BVP.lsNIPCE.ST.mean_uy_active_pc;
mean_u = BVP.lsNIPCE.ST.mean_u_pc;

% Hyperparameters and covariance matrices
rho_x = BVP.hyp.ST.rho(1); % rho for X dimension
rho_y = BVP.hyp.ST.rho(2); % rho for Y dimension
C_e = BVP.hyp.ST.C_e;
C_d_x = BVP.hyp.ST.C_d{1}; % Covariance for X dimension
C_d_y = BVP.hyp.ST.C_d{2}; % Covariance for Y dimension

% St. Venant material uncertainty
sigE = BVP.UQ.sig_E;

%% Initialize Posterior Mean and Covariance
mean_u_y = zeros(GDOFs, 1);
cov_u_y = zeros(GDOFs, GDOFs);

% Stabilize covariance matrices
cov_ux_red = cov_ux_active +1e-10 * sigE * eye(size(activeDOFsX, 2));
cov_uy_red = cov_uy_active +1e-10 * sigE * eye(size(activeDOFsY, 2));

%% Compute Posterior for Displacement in X
cov_ux_y_red = pinv(rho_x ^ 2 * n_0 * Px' * pinv(C_d_x + C_e) * Px + pinv(cov_ux_red));
mean_ux_y_red = cov_ux_y_red * (rho_x * Px' * pinv(C_d_x + C_e) * sum(BVP.obs.Y_exp_x(:, 1:n_0), 2) + pinv(cov_ux_red) * mean_ux_red);

% Scale posterior mean and covariance
mean_ux_y_red = rho_x * mean_ux_y_red;
cov_ux_y_red = rho_x ^ 2 * cov_ux_y_red;

% True response in X
mean_UZx_y_red = Px * mean_ux_y_red;
cov_UZx_y_red = Px * cov_ux_y_red * Px' + C_d_x;
ci_UZx_y_red = sqrt(diag(cov_UZx_y_red)) * 1.96;

%% Compute Posterior for Displacement in Y
cov_uy_y_red = pinv(rho_y ^ 2 * n_0 * Py' * pinv(C_d_y + C_e) * Py + pinv(cov_uy_red));
mean_uy_y_red = cov_uy_y_red * (rho_y * Py' * pinv(C_d_y + C_e) * sum(BVP.obs.Y_exp_y(:, 1:n_0), 2) + pinv(cov_uy_red) * mean_uy_red);

% Scale posterior mean and covariance
mean_uy_y_red = rho_y * mean_uy_y_red;
cov_uy_y_red = rho_y ^ 2 * cov_uy_y_red;

% True response in Y
mean_UZy_y_red = Py * mean_uy_y_red;
cov_UZy_y_red = Py * cov_uy_y_red * Py' + C_d_y;
ci_UZy_y_red = sqrt(diag(cov_UZy_y_red)) * 1.96;

%% Update Posterior Displacements
mean_u_y(activeDOFsX) = mean_ux_y_red; % X-direction from statFEM
mean_u_y(prescribedDOFsX) = mean_u(prescribedDOFsX); % X-direction from SFEM
mean_u_y(activeDOFsY) = mean_uy_y_red; % Y-direction from statFEM
mean_u_y(prescribedDOFsY) = mean_u(prescribedDOFsY); % Y-direction from SFEM
mean_u_y_bottomNodesDofx = mean_u_y(bottomNodesDofx);

cov_u_y(activeDOFsX, activeDOFsX) = cov_ux_y_red;
cov_u_y(activeDOFsY, activeDOFsY) = cov_uy_y_red;
cov_u_y = 0.5 * (cov_u_y + cov_u_y'); % Ensure symmetry
ci_cov_u_y_bottomNodesDofx = sqrt(diag(cov_u_y(bottomNodesDofx, bottomNodesDofx))) * 1.96;

% Surface computations
cov_u_y_x = cov_u_y(1:2:end, 1:2:end);
std_u_y_x = sqrt(diag(cov_u_y_x));
std_u_y_x_Surf = makeSurf(elementNodes, std_u_y_x);

mean_ux_y = mean_u_y(1:2:end);
mean_uy_y = mean_u_y(2:2:end);
mean_ux_y_Surf = makeSurf(elementNodes, mean_ux_y);

%% Assign Results Back to BVP
BVP.statFEM.ST.nAssimilation = n_0;
BVP.statFEM.ST.mean_u_y = mean_u_y;
BVP.statFEM.ST.mean_ux_y = mean_ux_y;
BVP.statFEM.ST.mean_uy_y = mean_uy_y;
BVP.statFEM.ST.mean_ux_y_Surf = mean_ux_y_Surf;
BVP.statFEM.ST.cov_u_y = cov_u_y;
BVP.statFEM.ST.std_u_y_x_Surf = std_u_y_x_Surf;
BVP.statFEM.ST.mean_u_y_bottomNodesDofx = mean_u_y_bottomNodesDofx;
BVP.statFEM.ST.ci_cov_u_y_bottomNodesDofx = ci_cov_u_y_bottomNodesDofx;

BVP.statFEM.ST.mean_UZx_y_red = mean_UZx_y_red;
BVP.statFEM.ST.ci_UZx_y_red = ci_UZx_y_red;
BVP.statFEM.ST.mean_UZy_y_red = mean_UZy_y_red;
BVP.statFEM.ST.ci_UZy_y_red = ci_UZy_y_red;
BVP.statFEM.ST.n_0 = n_0;

end
