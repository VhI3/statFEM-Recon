function BVP = infinitePlate_FEM_processLE_monteCarlo(BVP)
% INFINITEPLATE_FEM_PROCESSLE_MONTECARLO
% Perform Monte Carlo simulation for the Linear Elastic model
%
% Inputs:
%   BVP - Boundary value problem structure containing mesh, material, and
%         boundary condition data.
%
% Outputs:
%   BVP - Updated structure with Monte Carlo simulation results.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
activeDOFs = BVP.BC.activeDOFs;
prescribedDOFs = BVP.BC.prescribedDOFs;
force = BVP.BC.force;
nMCS = BVP.UQ.nMCS;
GDOFs = BVP.msh.GDOFs;
E_MC = BVP.mc.E_MC; % Monte Carlo samples for Young's modulus
elementNodes = BVP.msh.elementNodes;
bottomNodesDofx = BVP.msh.bottomNodesDofx;

% Preallocate displacement array
u_mc = zeros(nMCS, GDOFs);

% Progress bar
pb = CmdLineProgressBar('Processing MC for Linear Elastic...');

%% Monte Carlo Simulation
for ii = 1:nMCS
  pb.print(ii, nMCS);
  
  % Update material properties
  BVP.material.E = E_MC(ii);
  
  % Calculate global stiffness matrix
  stiffness = globalstiffness_linearElastic2D(BVP);
  
  % Partition stiffness matrix
  Krr = sparse(stiffness(activeDOFs, activeDOFs));
  Kru = sparse(stiffness(activeDOFs, prescribedDOFs));
  Kur = sparse(stiffness(prescribedDOFs, activeDOFs));
  Kuu = sparse(stiffness(prescribedDOFs, prescribedDOFs));
  
  % Static analysis
  Rr = sparse(force(activeDOFs)); % Reduced force vector
  Uu = sparse(zeros(length(prescribedDOFs), 1)); % Prescribed displacement
  Ur = Krr \ (Rr - Kru * Uu); % Solve for active DOFs
  Ru = Kuu * Uu + Kur * Ur; % Reactions at prescribed DOFs
  
  % Store displacement
  u = zeros(GDOFs, 1);
  u(activeDOFs) = Ur;
  u(prescribedDOFs) = Uu;
  u_mc(ii, :) = u;
end

%% Statistical Analysis
mean_u_mc = mean(u_mc, 1); % Mean displacement
std_u_mc = std(u_mc, 0, 1); % Standard deviation of displacement

% Separate x and y components
mean_ux_mc = mean_u_mc(1:2:end);
mean_uy_mc = mean_u_mc(2:2:end);
std_ux_mc = std_u_mc(1:2:end);
std_uy_mc = std_u_mc(2:2:end);

% Generate surface representations
mean_ux_mc_Surf = makeSurf(elementNodes, mean_ux_mc);
mean_uy_mc_Surf = makeSurf(elementNodes, mean_uy_mc);
std_ux_mc_Surf = makeSurf(elementNodes, std_ux_mc);
std_uy_mc_Surf = makeSurf(elementNodes, std_uy_mc);

% Displacement at bottom nodes
mean_ux_bottomNode_mc = mean_u_mc(bottomNodesDofx);
std_ux_bottomNode_mc = std_u_mc(bottomNodesDofx);

%% Assign Back to BVP
BVP.mc.LE.u_mc = u_mc;
BVP.mc.LE.mean_u_mc = mean_u_mc;
BVP.mc.LE.std_u_mc = std_u_mc;
BVP.mc.LE.mean_ux_mc = mean_ux_mc;
BVP.mc.LE.mean_uy_mc = mean_uy_mc;
BVP.mc.LE.std_ux_mc = std_ux_mc;
BVP.mc.LE.std_uy_mc = std_uy_mc;
BVP.mc.LE.mean_ux_mc_Surf = mean_ux_mc_Surf;
BVP.mc.LE.mean_uy_mc_Surf = mean_uy_mc_Surf;
BVP.mc.LE.std_ux_mc_Surf = std_ux_mc_Surf;
BVP.mc.LE.std_uy_mc_Surf = std_uy_mc_Surf;
BVP.mc.LE.mean_ux_bottomNode_mc = mean_ux_bottomNode_mc;
BVP.mc.LE.std_ux_bottomNode_mc = std_ux_bottomNode_mc;

end
