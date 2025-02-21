function BVP = infinitePlate_balancedStatus_ST(BVP)
% INFINITEPLATE_BALANCEDSTATUS_ST Computes the balanced equation status for nonlinear FEM analysis.
%   This function evaluates the balance of forces for a finite element
%   model in a nonlinear steady-state analysis and interpolates results
%   for visualization.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       Mesh data (element connectivity) and FEM results (balanced equations).
%
% Outputs:
%   BVP - Updated structure with computed balanced equation components and norms.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
elementNodes = BVP.msh.elementNodes; % Element connectivity
balancedEq = BVP.proc.ST.b; % Balanced equation vector

%% Calculation
% Extract balanced forces in x and y directions
balancedEq_x = balancedEq(1:2:end); % X-component of balanced forces
balancedEq_y = balancedEq(2:2:end); % Y-component of balanced forces

% Interpolate balanced forces to the surface for visualization
balancedEq_x_Surf = makeSurf(elementNodes, balancedEq_x); % Surface interpolation for X
balancedEq_y_Surf = makeSurf(elementNodes, balancedEq_y); % Surface interpolation for Y

% Compute the norm of the balanced equation vector
norm_balancedEq = norm(balancedEq);

%% Assign Back to BVP
BVP.proc.ST.balancedEq = balancedEq; % Full balanced equation vector
BVP.proc.ST.balancedEq_x = balancedEq_x; % Balanced forces in X-direction
BVP.proc.ST.balancedEq_y = balancedEq_y; % Balanced forces in Y-direction
BVP.proc.ST.balancedEq_x_Surf = balancedEq_x_Surf; % Surface interpolation for X
BVP.proc.ST.balancedEq_y_Surf = balancedEq_y_Surf; % Surface interpolation for Y
BVP.proc.ST.norm_balancedEq = norm_balancedEq; % Norm of the balanced equation vector

end
