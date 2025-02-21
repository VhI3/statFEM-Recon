function BVP = infinitePlate_balancedStatus_LE(BVP)
% INFINITEPLATE_BALANCEDSTATUS_LE Computes the balanced equation status for a linear elastic plate.
%   This function calculates the balance of forces to verify equilibrium in
%   a finite element model for a linear elastic infinite plate.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       msh.GDOFs              - Global degrees of freedom
%       msh.elementNodes       - Element connectivity
%       BC.activeDOFs          - Active degrees of freedom
%       BC.prescribedDOFs      - Prescribed degrees of freedom
%       proc.LE.Krr, Kru, Kur, Kuu - Stiffness matrix partitions
%       proc.LE.Rr, Ru         - Force vectors
%       proc.LE.Ur, Uu         - Displacement vectors
%
% Outputs:
%   BVP - Updated BVP structure with the balanced equation results.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
GDOFs = BVP.msh.GDOFs; % Total global degrees of freedom
elementNodes = BVP.msh.elementNodes; % Element connectivity
activeDOFs = BVP.BC.activeDOFs; % Active degrees of freedom
prescribedDOFs = BVP.BC.prescribedDOFs; % Prescribed degrees of freedom

% Stiffness matrix partitions
Krr = BVP.proc.LE.Krr;
Kru = BVP.proc.LE.Kru;
Kur = BVP.proc.LE.Kur;
Kuu = BVP.proc.LE.Kuu;

% Force and displacement vectors
Rr = BVP.proc.LE.Rr;
Ru = BVP.proc.LE.Ru;
Uu = BVP.proc.LE.Uu;
Ur = BVP.proc.LE.Ur;

%% Calculation
% Initialize the balanced equation vector
balancedEq = zeros(GDOFs, 1);

% Compute balance of forces for active and prescribed DOFs
balancedEq(activeDOFs) = Krr * Ur + Kru * Uu - Rr; % Active DOFs
balancedEq(prescribedDOFs) = Kur * Ur + Kuu * Uu - Ru; % Prescribed DOFs

% Separate balanced forces into x and y components
balancedEq_x = balancedEq(1:2:end); % X-direction
balancedEq_y = balancedEq(2:2:end); % Y-direction

% Interpolate balanced forces to the surface for visualization
balancedEq_x_Surf = makeSurf(elementNodes, balancedEq_x); % Surface in X
balancedEq_y_Surf = makeSurf(elementNodes, balancedEq_y); % Surface in Y

% Compute the norm of the balanced equation vector
norm_balancedEq = norm(balancedEq);

%% Assign back to BVP
BVP.proc.LE.balancedEq = balancedEq; % Full balanced equation vector
BVP.proc.LE.balancedEq_x = balancedEq_x; % Balanced forces in X-direction
BVP.proc.LE.balancedEq_y = balancedEq_y; % Balanced forces in Y-direction
BVP.proc.LE.balancedEq_x_Surf = balancedEq_x_Surf; % Surface interpolation in X
BVP.proc.LE.balancedEq_y_Surf = balancedEq_y_Surf; % Surface interpolation in Y
BVP.proc.LE.norm_balancedEq = norm_balancedEq; % Norm of the balanced equation vector

end
