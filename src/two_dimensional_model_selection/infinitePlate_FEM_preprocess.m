function BVP = infinitePlate_FEM_preprocess(BVP)
% INFINITEPLATE_FEM_PREPROCESS Preprocesses FEM data for a plate with a hole.
%
%   BVP = INFINITEPLATE_FEM_PREPROCESS(BVP) imports mesh data, calculates geometry,
%   material properties, boundary conditions, and prepares data for finite element
%   analysis, uncertainty quantification, and observation handling.
%
% Inputs:
%   BVP - Structure containing necessary inputs and fields for FEM preprocessing.
%
% Outputs:
%   BVP - Updated structure with calculated properties and settings.
%
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)
% Project: statFEM-Recon

%% Import and Assign Mesh Data
BVP = importMesh(BVP); % Import mesh
nodeCoordinates = BVP.msh.nodeCoordinates; % Node coordinates
elementNodes = BVP.msh.elementNodes; % Element connectivity

%% Geometry Properties
L = 4; H = 4; % Plate dimensions (length, height) [m]
R = nodeCoordinates(1, 2); % Radius of the hole [m]
T = 1; % Thickness of the plate [m]
DIM = 2; % Problem dimension

%% Mesh Properties
nuElm = size(elementNodes, 1); % Number of elements
nuNodes = size(nodeCoordinates, 1); % Number of nodes
DOFs = 2; % Degrees of freedom per node
nElNode = 4; % Nodes per element
eDOFs = DOFs * nElNode; % Degrees of freedom per element
GDOFs = nuNodes * DOFs; % Global degrees of freedom

% Generate element DOFs
elementDOFs = zeros(nuElm, eDOFs);
elementDOFs(:, 1:2:eDOFs) = 2 * elementNodes - 1;
elementDOFs(:, 2:2:eDOFs) = 2 * elementNodes;

% Node sets for boundary conditions
leftNodes = find(abs(nodeCoordinates(:, 1)) < eps);
rightNodes = find(abs(nodeCoordinates(:, 1) - L) < eps);
topNodes = find(abs(nodeCoordinates(:, 2) - H) < eps);
bottomNodes = find(abs(nodeCoordinates(:, 2)) < eps);

% Corner nodes
cornerNode = intersect(rightNodes, bottomNodes);
cornerNodeDofX = DOFs * cornerNode - 1;

% Bottom nodes ordered for displacement
bottomNodesDofx = 2 * bottomNodes - 1;
bottomNodesCoor = nodeCoordinates(bottomNodes, 1);

%% Material Properties
E = 200; % Young's modulus
nu = 0.25; % Poisson's ratio
lamda = @(E) nu * E / (1 + nu) / (1 - 2 * nu); % Lame's constant 1
mu = @(E) E / (2 * (1 + nu)); % Lame's constant 2
Dmatrix = @(E) [lamda(E) + 2 * mu(E), lamda(E), 0;
  lamda(E), lamda(E) + 2 * mu(E), 0;
  0, 0, mu(E)];

% Constitutive matrices
C_pstress = 1 / (1 - nu ^ 2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
C_pstrain = 1 / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, 1/2 - nu];

%% Uncertainty Quantification
mu_E = 200; % Mean Young's modulus
sig_E = 0.1 * mu_E; % Standard deviation of Young's modulus
nMCS = 2; % Monte Carlo samples

%% Processing Properties
u = zeros(GDOFs, 1); % Initialize displacement vector
BVP.proc.LE.u = zeros(GDOFs, 1); % Linear elasticity initial guess
xiRange1 = [-1, 1]; % Gauss integration range (can be modified)
xiRange2 = [-1, 1];
NGT = 2; % Number of Gauss points
[xi1, w1] = Gauss_int(NGT); % Gauss points and weights for xi1
[xi2, w2] = Gauss_int(NGT); % Gauss points and weights for xi2
dispSF = 1; % Displacement scaling factor for visualization

%% Boundary Conditions
% Dirichlet boundary conditions
prescribedDOFsX = DOFs * leftNodes - 1; % X-direction constraints
prescribedDOFsY = DOFs * bottomNodes; % Y-direction constraints
prescribedDOFs = sort([prescribedDOFsX; prescribedDOFsY]);

% Prescribed displacements and forces
dbc = u(prescribedDOFs); % Displacement boundary conditions
ndbc = size(prescribedDOFs, 1); % Number of constrained DOFs
traction = 100; % Applied traction

% Correct right nodes for ordering
if size(rightNodes, 1) > 2
  rightNodesCorr = [rightNodes(1); rightNodes(3:end); rightNodes(2)];
else
  rightNodesCorr = rightNodes;
end

% Right edge for force application
rightEdge_y1 = nodeCoordinates(rightNodesCorr, 2);
rightEdge_y = [rightEdge_y1(1:end - 1), rightEdge_y1(2:end)];
rightEdge_y_dofX1 = DOFs * rightNodesCorr - 1;
rightEdge_y_dofX = [rightEdge_y_dofX1(1:end - 1), rightEdge_y_dofX1(2:end)];
force = zeros(GDOFs, 1);

for ey = 1:size(rightEdge_y, 1)
  ed = rightEdge_y_dofX(ey, :);
  ec = rightEdge_y(ey, :);
  force(ed, 1) = force(ed, 1) + abs(ec(2) - ec(1)) / 2 * traction * ones(2, 1);
end

tbc = [];
tbc = [tbc rightNodes];
ntbc = size(tbc, 1);

% Active and inactive DOFs
activeDOFsX = setdiff(1:2:GDOFs, prescribedDOFsX);
activeDOFsY = setdiff(2:2:GDOFs, prescribedDOFsY);
activeDOFs = setdiff(1:GDOFs, prescribedDOFs);

%% Newton-Raphson Method
tol = 1e-10; % Convergence tolerance
maxit = 20; % Maximum iterations
reit = 0; % Reduction index for step size
maxreit = 6; % Maximum reduction steps
timeInterval = 1; % Initial time interval

%% Non-Intrusive Least Square PCE
p_order = 8; % Polynomial order for PCE
[~, Psi_s, ~, PsiSqNorm, P_PCE] = Hermite_PC(1, p_order); % Hermite basis
N_s = 2 * P_PCE; % Sample size for PCE

% Log-normal transformation for uncertainty quantification
lambda = log(mu_E ^ 2 / sqrt(mu_E ^ 2 + sig_E ^ 2));
zeta = sqrt(log(1 + sig_E ^ 2 / mu_E ^ 2));
rng(3); uniform_PC = rand(N_s, 1); % Uniform samples for PCE
rng(3); uniform_MC = rand(nMCS, 1); % Uniform samples for MC

Xi_PC = norminv(uniform_PC, 0, 1); % Gaussian samples for PCE
Xi_MC = norminv(uniform_MC, 0, 1); % Gaussian samples for MC

% Compute Young's modulus samples
E_PC = zeros(N_s, 1);
E_MC = zeros(nMCS, 1);
P_E = 3;

for i = 1:P_E
  
  Xi_PC_rational = arrayfun(@(x) sym(sprintf('%.17g', x)), Xi_PC(:, 1));
  bigPsi_PC = (1 / sqrt(factorial(i))) * double(subs(Psi_s{i, 1}, Xi_PC_rational));
  %bigPsi_PC = (1 / sqrt(factorial(i))) * double(subs(Psi_s{i, 1}, Xi_PC(:, 1)));
  E_PC = E_PC + mu_E * (zeta ^ (i - 1) / factorial(i - 1)) * bigPsi_PC;
  
  Xi_MC_rational = arrayfun(@(x) sym(sprintf('%.17g', x)), Xi_MC(:, 1));
  bigPsi_MC = (1 / sqrt(factorial(i))) * double(subs(Psi_s{i, 1}, Xi_MC_rational));
  %bigPsi_MC = (1 / sqrt(factorial(i))) * double(subs(Psi_s{i, 1}, Xi_MC(:, 1)));
  E_MC = E_MC + mu_E * (zeta ^ (i - 1) / factorial(i - 1)) * bigPsi_MC;
end

% -------------------- 10. Observation data -------------------------------
% sensorLocation =nodeCoordinates; save 'sensorLocation.mat' sensorLocation
load 'sensorLocation.mat' sensorLocation
[sensorIndex, ~] = dsearchn(nodeCoordinates, sensorLocation);
senCoor = [nodeCoordinates(sensorIndex, 1) nodeCoordinates(sensorIndex, 2)];
bottomSen = find(abs(senCoor(:, 2)) - eps < 0);
[~, BB] = sort(senCoor(bottomSen, 1));
bottomSen = bottomSen(BB);
bottomSen_DofX = 2 * bottomSen - 1;
% hold on
% scatter(senCoor(:,1),senCoor(:,2),'SizeData',4,'MarkerFaceColor','k','MarkerEdgeColor','k',...
%     'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',1);
% hold on
% plot([0.4;4],[0; 0],'r','LineWidth',1.5);
% % % hold on
% % % scatter(senCoor(bottomSen,1),senCoor(bottomSen,2),'r*');

%% Case 1
nSen = size(senCoor, 1);
nSam = 100;
epsi = 4e-4;
l_d_sample = 2;
sig_d_sample = 0.2;
rho_sample = 1.5;
C_e_sample = epsi * eye(nSen, nSen);
e_sample = mvnrnd(zeros(nSen, 1), C_e_sample, nSam)';
C_d_sample = sqexp(senCoor, senCoor, log(sig_d_sample), log(l_d_sample));
d_sample = mvnrnd(zeros(nSen, 1), C_d_sample, nSam)';
%
%% Assign back to BVP
% -------------------- 1. Geometry properties -----------------------------
BVP.geometry.L = L;
BVP.geometry.H = H;
BVP.geometry.T = T;
BVP.geometry.R = R;
BVP.geometry.DIM = DIM;
% -------------------- 2. Mesh properties ---------------------------------
BVP.msh.nElm = nuElm;
BVP.msh.nodeCoordinates = nodeCoordinates;
BVP.msh.DOFs = DOFs;
BVP.msh.elementNodes = elementNodes;
BVP.msh.elementDOFs = elementDOFs;
BVP.msh.nuNodes = nuNodes;
BVP.msh.eDOFs = eDOFs;
BVP.msh.GDOFs = GDOFs;
BVP.msh.leftNodes = leftNodes;
BVP.msh.rightNodes = rightNodes;
BVP.msh.topNodes = topNodes;
BVP.msh.bottomNodes = bottomNodes;
BVP.msh.bottomNodesDofx = bottomNodesDofx;
BVP.msh.bottomNodesCoor = bottomNodesCoor;
BVP.msh.cornerNode = cornerNode;
BVP.msh.cornerNodeDofX = cornerNodeDofX;
% -------------------- 3. Material properties -----------------------------
BVP.material.E = E;
BVP.material.nu = nu;
BVP.material.lamda = lamda;
BVP.material.mu = mu;
BVP.material.Dmatrix = Dmatrix;
BVP.material.C_pstrain = C_pstrain;
BVP.material.C_pstress = C_pstress;
% -------------------- 4. Uncertain material properties -------------------
BVP.UQ.mu_E = mu_E;
BVP.UQ.sig_E = sig_E;
BVP.UQ.nMCS = nMCS;
% -------------------- 6. processing properties ---------------------------
% BVP.proc.u                   = u;
BVP.proc.xiRange1 = xiRange1;
BVP.proc.xiRange2 = xiRange2;
BVP.proc.dispSF = dispSF;
BVP.proc.GP.NGT = NGT;
BVP.proc.GP.xi1 = xi1;
BVP.proc.GP.xi2 = xi2;
BVP.proc.GP.w11 = w1;
BVP.proc.GP.w12 = w2;
% -------------------- 7. Boundary conditions -----------------------------
BVP.BC.prescribedDOFs = prescribedDOFs;
BVP.BC.prescribedDOFsX = prescribedDOFsX;
BVP.BC.prescribedDOFsY = prescribedDOFsY;
BVP.BC.activeDOFs = activeDOFs;
BVP.BC.activeDOFsX = activeDOFsX;
BVP.BC.activeDOFsY = activeDOFsY;
BVP.BC.ndbc = ndbc;
BVP.BC.ntbc = ntbc;
BVP.BC.dbc = dbc;
BVP.BC.force = force;
% -------------------- 8. newton raphson ----------------------------------
BVP.nRaphson.tol = tol;
BVP.nRaphson.maxit = maxit;
BVP.nRaphson.reit = reit;
BVP.nRaphson.maxreit = maxreit;
BVP.nRaphson.timeInterval = timeInterval;
% -------------------- 9. non-intrusive least square PCE ------------------
BVP.lsNIPCE.N_s = N_s;
BVP.lsNIPCE.P_PCE = P_PCE;
BVP.lsNIPCE.Psi_s = Psi_s;
BVP.lsNIPCE.PsiSqNorm = PsiSqNorm;
BVP.lsNIPCE.zeta = zeta;
BVP.lsNIPCE.lambda = lambda;
BVP.lsNIPCE.uniform_PC = uniform_PC;
BVP.mc.uniform_MC = uniform_MC;
BVP.lsNIPCE.Xi_PC = Xi_PC;
BVP.mc.Xi_MC = Xi_MC;
BVP.lsNIPCE.P_E = P_E;
BVP.lsNIPCE.E_PC = E_PC;
BVP.mc.E_MC = E_MC;
% -------------------- 10. Observation data -------------------------------
BVP.obs.senCoor = senCoor;
BVP.obs.nSam = nSam;
BVP.obs.nSen = nSen;
BVP.obs.l_d_sample = l_d_sample;
BVP.obs.sig_d_sample = sig_d_sample;
BVP.obs.rho_sample = rho_sample;
BVP.obs.epsi = epsi;
BVP.obs.C_d_sample = C_d_sample;
BVP.obs.d_sample = d_sample;
BVP.obs.e_sample = e_sample;
BVP.obs.C_e_sample = C_e_sample;
BVP.obs.sensorIndex = sensorIndex;
BVP.obs.bottomSen_DofX = bottomSen_DofX;
BVP.obs.bottomSen = bottomSen;
end
