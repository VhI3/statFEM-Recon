function BVP = tensionBar_1D_FEM_processLE(BVP)
% tensionBar_1D_FEM_processLE
% Finite Element Method analysis for 1D tension bar using linear elasticity.
% This function computes the displacement based on
% material properties and external loading.
% --------------------------------------------------------------
% Inputs:
%   BVP - Boundary Value Problem structure containing material properties,
%         geometry, and loading.
% Outputs:
%   BVP - Updated BVP structure with FEM results.
% --------------------------------------------------------------
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
f_bar = BVP.loading.f_bar; % Applied tip load
A = BVP.geometry.A; % Cross-sectional area
nElm = BVP.mesh.nElm; % Number of elements
nodeCoordinates = BVP.mesh.nodeCoordinates; % Node coordinates
E = BVP.material.mu_E;
activeDOFs = BVP.fem.activeDOFs;
GDof = BVP.fem.GDof;
elementNodes = BVP.fem.elementNodes;
%% FEM Processing
displacement = FEM_Bar_deter_Tipload(E, f_bar, A, nElm, elementNodes, nodeCoordinates, activeDOFs, GDof);
%% Assign results back to BVP
BVP.results.LE.displacement = displacement;

disp('2. FEM Linear Elasticity Analysis Completed')
end
