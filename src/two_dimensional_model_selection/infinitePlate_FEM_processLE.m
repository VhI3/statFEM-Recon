function BVP = infinitePlate_FEM_processLE(BVP)
% INFINITEPLATEFEMPROCESSLE Performs FEM processing for a linear elastic infinite plate.
%   This function calculates the displacement and force vectors for a
%   2D linear elastic problem using the finite element method (FEM).
%
% Inputs:
%   BVP - Boundary value problem structure containing mesh, boundary
%         conditions, and processing fields.
%
% Outputs:
%   BVP - Updated structure with computed displacement, force, and
%         stiffness matrix components.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
activeDOFs = BVP.BC.activeDOFs; % Active degrees of freedom
prescribedDOFs = BVP.BC.prescribedDOFs; % Prescribed degrees of freedom
force = BVP.BC.force; % Force vector
u = BVP.proc.LE.u; % Displacement vector
elementNodes = BVP.msh.elementNodes; % Element connectivity
bottomNodesDofx = BVP.msh.bottomNodesDofx; % Bottom nodes DOFs (x-direction)

%% Calculation of the system stiffness matrix
% Compute the global stiffness matrix for the 2D linear elastic problem
stiffness = globalstiffness_linearElastic2D(BVP);

% Partition the stiffness matrix
Krr = sparse(stiffness(activeDOFs, activeDOFs)); % Active DOFs stiffness
Kru = sparse(stiffness(activeDOFs, prescribedDOFs)); % Active-Prescribed coupling
Kur = sparse(stiffness(prescribedDOFs, activeDOFs)); % Prescribed-Active coupling
Kuu = sparse(stiffness(prescribedDOFs, prescribedDOFs)); % Prescribed DOFs stiffness

%% Solve the system
% Static analysis
Rr = sparse(force(activeDOFs)); % Reduced force vector
Uu = sparse(u(prescribedDOFs)); % Prescribed displacements
Ur = Krr \ (Rr - Kru * Uu); % Solve for unknown displacements
Ru = Kuu * Uu + Kur * Ur; % Compute reaction forces

% Update displacement and force vectors
u(activeDOFs) = Ur;
force(prescribedDOFs) = Ru;

%% Post-processing
% Separate displacements into x and y components
ux = u(1:2:end, :); % Displacement in x-direction
uy = u(2:2:end, :); % Displacement in y-direction

% Interpolate displacements to the surface
ux_Surf = makeSurf(elementNodes, ux); % Surface displacement in x
uy_Surf = makeSurf(elementNodes, uy); % Surface displacement in y

% Bottom node displacement in x-direction
ux_bottomNode = u(bottomNodesDofx);

%% Assign back to BVP
BVP.proc.LE.u = u;
BVP.proc.LE.ux = ux;
BVP.proc.LE.uy = uy;
BVP.proc.LE.ux_Surf = ux_Surf;
BVP.proc.LE.uy_Surf = uy_Surf;
BVP.proc.LE.ux_bottomNode = ux_bottomNode;
BVP.proc.LE.Krr = Krr;
BVP.proc.LE.Kru = Kru;
BVP.proc.LE.Kur = Kur;
BVP.proc.LE.Kuu = Kuu;
BVP.proc.LE.Rr = Rr;
BVP.proc.LE.Ru = Ru;
BVP.proc.LE.Uu = Uu;
BVP.proc.LE.Ur = Ur;
BVP.proc.LE.force = force;
end
