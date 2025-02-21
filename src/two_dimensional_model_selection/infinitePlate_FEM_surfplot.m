function BVP = infinitePlate_FEM_surfplot(BVP, flag)
%
% infinitePlate_FEM_surfplot - Generates surface plots for various displacement and stress fields.
%
% This function visualizes the results of different analysis methods applied to the infinite plate problem.
% It includes:
%   - Monte Carlo (MC) simulation results
%   - Polynomial Chaos (PC) expansion results
%   - Statistical Finite Element Method (statFEM) posterior results
%   - Residual and stress fields for both Linear Elastic (LE) and St. Venant Kirchhoff (ST) models
%
% Inputs:
%   BVP  - Boundary Value Problem structure containing all simulation results.
%   flag - Plot control flag (1 to generate plots, 0 to skip plotting).
%
% Outputs:
%   BVP  - Updated BVP structure (if any modifications are added in future).
%
% Usage:
%   BVP = infinitePlate_FEM_surfplot(BVP, 1); % Generates and displays plots
%
% Note:
% - The function uses subplots to compare different methods side-by-side.
% - Compatible with both MATLAB and GNU Octave.
% - Ensure 'BVP' structure has valid computed fields before plotting.
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Assign from BVP
%%
if flag == 1
    close all;
    %%
    subplot(3, 8, 1)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.mc.LE.mean_ux_mc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
        colorbar('eastoutside'); % Octave-compatible colorbar
    else
        colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x^{MC}$ LE')
    %%
    subplot(3, 8, 2)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.lsNIPCE.LE.mean_ux_pc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
        colorbar('eastoutside'); % Octave-compatible colorbar
    else
        colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x^{PC}$ LE')
    %%
    subplot(3, 8, 3)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.mc.LE.std_ux_mc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
        colorbar('eastoutside'); % Octave-compatible colorbar
    else
        colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x^{mC}$ LE')
    %%
    subplot(3, 8, 4)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.lsNIPCE.LE.std_ux_pc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
        colorbar('eastoutside'); % Octave-compatible colorbar
    else
        colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x^{pC}$ LE')
    %%
    subplot(3, 8, 5)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.postproc.LE.sigma_xx_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('$\sigma_{xx}$ LE')
    %%
    subplot(3, 8, 6)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.LE.mean_ux_y_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x|y$ LE')
    %%
    subplot(3, 8, 7)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.LE.std_u_y_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x|y$ LE')
    %%
    subplot(3, 8, 8)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.postproc.LE.sigma_xx_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('$\sigma_{xx}| Y$ LE')
    %%
    subplot(3, 8, 9)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.mc.ST.mean_ux_mc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x^{MC}$ ST')
    %%
    subplot(3, 8, 10)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.lsNIPCE.ST.mean_ux_pc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x^{PC}$ ST')
    %%
    subplot(3, 8, 11)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.mc.ST.std_ux_mc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x^{mC}$ ST')
    %%
    subplot(3, 8, 12)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.lsNIPCE.ST.std_ux_pc_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x^{pC}$ ST')
    %%
    subplot(3, 8, 13)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.postproc.ST.pk1Stress_xx_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('$1PK_{xx}$ ST')
    %%
    subplot(3, 8, 14)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.ST.mean_ux_y_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('mean $u_x|y$ ST')
    %%
    subplot(3, 8, 15)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.ST.std_u_y_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('std $u_x|y$ ST')
    %%
    subplot(3, 8, 16)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.postproc.ST.pk1Stress_xx_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('$\sigma_{xx}| Y$ ST')
    %%
    subplot(3, 8, 17)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.proc.LE.balancedEq_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('Prior Residuum in $x$-direction LE')
    %%
    subplot(3, 8, 18)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.LE.balancedEq_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('Posterior Residuum in $x$-direction LE')
    %%
    subplot(3, 8, 19)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.proc.ST.balancedEq_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('Prior Residuum in $x$-direction SV')
    %%
    subplot(3, 8, 20)
    patch(BVP.msh.xSurf, BVP.msh.ySurf, BVP.statFEM.ST.balancedEq_x_Surf, 'FaceColor', 'interp');
    axis equal; view(0, 90); colormap(jet);
    
    if exist('OCTAVE_VERSION', 'builtin') % Check if running in Octave
      colorbar('eastoutside'); % Octave-compatible colorbar
    else
      colorbar('vert'); % MATLAB-compatible colorbar
    end
    
    axis off; shading interp;
    title('Posterior Residuum in $x$-direction SV')
    
end

end
