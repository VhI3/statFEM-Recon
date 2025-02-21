function BVP = infinitePlate_FEM_processST_monteCarlo(BVP)
% INFINITEPLATE_FEM_PROCESSST_MONTECARLO
% Perform Monte Carlo simulation for the St. Venant Kirchhoff model.
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
nMCS = BVP.UQ.nMCS; % Number of Monte Carlo simulations
elementNodes = BVP.msh.elementNodes; % Element connectivity
E_MC = BVP.mc.E_MC; % Monte Carlo samples for Young's modulus
bottomNodesDofx = BVP.msh.bottomNodesDofx; % Bottom nodes for DOF x

% Preallocate displacement array
GDOFs = length(BVP.msh.nodeCoordinates) * 2; % Total DOFs
u_mc = zeros(nMCS, GDOFs);

% Initialize temporary BVP structure for modifications
BVP_tmp = BVP;

% Progress bar
pb = CmdLineProgressBar('MC of ST Processing...');

%% Monte Carlo Simulation
for i = 1:nMCS
    pb.print(i, nMCS);
    
    % Update material properties
    BVP_tmp.material.E = E_MC(i);
    
    % Compute displacements for the current realization
    BVP_tmp = infinitePlate_FEM_processST(BVP_tmp);
    u_mc(i, :) = BVP_tmp.proc.ST.u; % Store displacement vector
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
BVP.mc.ST.u_mc = u_mc;
BVP.mc.ST.mean_u_mc = mean_u_mc;
BVP.mc.ST.std_u_mc = std_u_mc;
BVP.mc.ST.mean_ux_mc = mean_ux_mc;
BVP.mc.ST.mean_uy_mc = mean_uy_mc;
BVP.mc.ST.std_ux_mc = std_ux_mc;
BVP.mc.ST.std_uy_mc = std_uy_mc;
BVP.mc.ST.mean_ux_mc_Surf = mean_ux_mc_Surf;
BVP.mc.ST.mean_uy_mc_Surf = mean_uy_mc_Surf;
BVP.mc.ST.std_ux_mc_Surf = std_ux_mc_Surf;
BVP.mc.ST.std_uy_mc_Surf = std_uy_mc_Surf;
BVP.mc.ST.mean_ux_bottomNode_mc = mean_ux_bottomNode_mc;
BVP.mc.ST.std_ux_bottomNode_mc = std_ux_bottomNode_mc;

end
