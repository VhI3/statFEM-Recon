function BVP = infinitePlate_FEM_processST(BVP)
% INFINITEPLATE_FEM_PROCESSST Performs FEM processing for nonlinear analysis of a plate.
%   This function implements a nonlinear Newton-Raphson iterative solver
%   for the finite element method (FEM) analysis of an infinite plate.
%
% Inputs:
%   BVP - Boundary value problem structure containing:
%       Mesh, boundary conditions, material properties, and FEM parameters.
%
% Outputs:
%   BVP - Updated structure containing computed displacements and related quantities.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
% General problem parameters
GDOFs = BVP.msh.GDOFs; % Total degrees of freedom
tol = BVP.nRaphson.tol; % Convergence tolerance
maxit = BVP.nRaphson.maxit; % Maximum iterations
reit = BVP.nRaphson.reit; % Reduction index
maxreit = BVP.nRaphson.maxreit; % Maximum refinement iterations
timeInterval = BVP.nRaphson.timeInterval; % Time step size

% Boundary conditions
ndbc = BVP.BC.ndbc; % Number of prescribed displacements
ntbc = BVP.BC.ntbc; % Number of traction boundary conditions
dbc = BVP.BC.dbc; % Displacement boundary conditions
prescribedDOFs = BVP.BC.prescribedDOFs; % Prescribed degrees of freedom
activeDOFs = BVP.BC.activeDOFs; % Active degrees of freedom
force = BVP.BC.force; % Force vector
bottomNodesDofx = BVP.msh.bottomNodesDofx; % Bottom node DOFs in x-direction

% Temporary structure for material and mesh properties
BVP_tmp.GP = BVP.proc.GP;
BVP_tmp.nElm = BVP.msh.nElm;
BVP_tmp.nodeCoordinates = BVP.msh.nodeCoordinates;
BVP_tmp.elementNodes = BVP.msh.elementNodes;
BVP_tmp.elementDOFs = BVP.msh.elementDOFs;
BVP_tmp.GDOFs = BVP.msh.GDOFs;
BVP_tmp.E = BVP.material.E; % Young's modulus
BVP_tmp.C = BVP.material.C_pstrain; % Plane strain matrix
BVP_tmp.nu = BVP.material.nu; % Poisson's ratio
BVP_tmp.lamda = BVP.material.lamda(BVP_tmp.E); % Lame parameter
BVP_tmp.mu = BVP.material.mu(BVP_tmp.E); % Shear modulus
BVP_tmp.Dmatrix = BVP.material.Dmatrix(BVP_tmp.E); % Elasticity matrix

%% Initialization
u = zeros(GDOFs, 1); % Nodal displacements
u_steps = u; % Displacement steps
cu = zeros(GDOFs, 1); % Converged nodal displacements
step = 0; % Step index
curtime = 0; % Current time
cnit = []; % Record iterative steps
counter = 0; % Iteration counter

%% Newton-Raphson Iteration
while curtime ~= 1
  counter = counter + 1;
  curtime = curtime + timeInterval;
  
  % Adjust time interval if exceeding final step
  if curtime > 1
    timeInterval = 1 - curtime + timeInterval;
    curtime = 1;
  end
  
  % Initialize error and iteration counters
  err = 1e6; % Initial error
  perr = err; % Previous error
  nit = 0; % Iteration count
  iu = zeros(GDOFs, 1); % Iterative displacements
  
  % Newton-Raphson loop
  while (err > tol) && (nit <= maxit)
    nit = nit + 1;
    
    % Compute stiffness matrix and residual force
    [k, r] = globalstiffness_stKirchhoff2D(BVP_tmp, u);
    
    % Apply boundary conditions
    f = zeros(GDOFs, 1); % External force vector
    
    if ntbc ~= 0
      f = force; % Apply traction boundary conditions
    end
    
    if ndbc ~= 0
      k(prescribedDOFs, :) = 0;
      k(prescribedDOFs, prescribedDOFs) = eye(ndbc);
      f(prescribedDOFs, :) = 0;
      
      if nit == 1
        f(prescribedDOFs, :) = dbc;
      end
      
    end
    
    % Compute right-hand side of governing equation
    b = curtime * f - r;
    
    if ndbc ~= 0
      b(prescribedDOFs) = curtime * dbc - u(prescribedDOFs);
    end
    
    % Solve for displacement increment
    du = k \ b;
    
    % Update displacement
    u = u + du;
    iu = iu + du;
    
    % Compute iterative error
    if nit > 1
      num = b(activeDOFs)' * b(activeDOFs);
      denom = 1 + f(activeDOFs)' * f(activeDOFs);
      err = num / denom;
    end
    
    % Handle divergence
    if err / perr > 1E6 && nit > 2
      nit = maxit + 1; % Mark as non-convergent
    else
      perr = err;
    end
    
  end
  
  % Check convergence
  if nit <= maxit
    % Converged
    reit = 0;
    step = step + 1;
    cu = u;
    cnit = [cnit, nit]; %#ok<AGROW>
    
    % Adjust time interval for faster convergence
    if length(cnit) >= 2 && all(cnit( end - 1:end) <= 5)
      timeInterval = timeInterval * 1.5;
    end
    
    % Save converged displacement
    u_steps = [u_steps, u]; %#ok<AGROW>
  else
    % Not converged
    if reit <= maxreit
      curtime = curtime - timeInterval; % Restore time step
      timeInterval = timeInterval / 4; % Refine time interval
      reit = reit + 1;
      u = cu; % Restore last converged displacement
    else
      return; % Stop analysis
    end
    
  end
  
end

%% Post-Processing
ux_steps = u_steps(1:2:end, :); % X-displacements over steps
uy_steps = u_steps(2:2:end, :); % Y-displacements over steps
u = u_steps(:, end); % Final displacements
ux = u(1:2:end); % Final X-displacement
uy = u(2:2:end); % Final Y-displacement
ux_Surf = makeSurf(BVP_tmp.elementNodes, ux); % Surface X-displacement
uy_Surf = makeSurf(BVP_tmp.elementNodes, uy); % Surface Y-displacement
ux_bottomNode = u(bottomNodesDofx); % Bottom node X-displacement

%% Assign Back to BVP
BVP.proc.ST.u_steps = u_steps;
BVP.proc.ST.ux_steps = ux_steps;
BVP.proc.ST.uy_steps = uy_steps;
BVP.proc.ST.u = u;
BVP.proc.ST.ux = ux;
BVP.proc.ST.uy = uy;
BVP.proc.ST.ux_Surf = ux_Surf;
BVP.proc.ST.uy_Surf = uy_Surf;
BVP.proc.ST.ux_bottomNode = ux_bottomNode;
BVP.proc.ST.b = b;

end
