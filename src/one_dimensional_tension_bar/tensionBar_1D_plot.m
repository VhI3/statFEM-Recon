function BVP = tensionBar_1D_plot(BVP)
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Assign Parameters
L = BVP.geometry.L;
nodeCoordinates = BVP.mesh.nodeCoordinates;
senCoor = BVP.obs.senCoor;
mu_u_pc = BVP.UQ.mu_u_pc;
C_u_pc = BVP.UQ.C_u_pc;
ci_u_pc = BVP.UQ.ci_u_pc;
u_mc = BVP.UQ.LE.u_mc;
pdf_u_tip_pc = BVP.UQ.LE.pdf_u_tip_pc;
cdf_u_tip_pc = BVP.UQ.LE.cdf_u_tip_pc;
xpc_pdf = BVP.UQ.LE.xpc_pdf;
xpc_cdf = BVP.UQ.LE.xpc_cdf;
pdf_u_tip_mc = BVP.UQ.LE.pdf_u_tip_mc;
cdf_u_tip_mc = BVP.UQ.LE.cdf_u_tip_mc;
xmc_pdf = BVP.UQ.LE.xmc_pdf;
xmc_cdf = BVP.UQ.LE.xmc_cdf;
mu_u_y = BVP.statFEM.mu_u_y;
mu_z = BVP.statFEM.mu_z;
C_u_y = BVP.statFEM.C_u_y;
C_z = BVP.statFEM.C_z;
ci_u_y = BVP.statFEM.ci_u_y;
ci_z = BVP.statFEM.ci_z;
nrep = BVP.obs.nrep;
y_obs = BVP.obs.y_obs;
close all;
%% figure 1: Prior problem
figure
plot_u_mean_pc = plot(nodeCoordinates, mu_u_pc, 'b-', 'LineWidth', 0.5);
set(plot_u_mean_pc, 'DisplayName', '$\mu_u$');
hold on
%
opts = {'EdgeColor', [0 0 1], 'FaceColor', [0 0 1], 'alpha', 0.5};
where = nodeCoordinates >- 1 & nodeCoordinates < nodeCoordinates(end) + 1;
plot_ci_c_u = fill_between(nodeCoordinates, (mu_u_pc - ci_u_pc), (mu_u_pc + ci_u_pc), where, opts{:});
set(plot_ci_c_u, 'DisplayName', '$95$\% CI')
grid on
%
hl = legend('show', 'FontSize', 8, 'location', 'northwest');
set(hl, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$u(x)$', 'Interpreter', 'latex')
axis square
ylim([-1 25]); xlim([-5 L + 5]);
title('Prior Dispplacement Response')
hold off
%
%% Subfigure 2: Monte Carlo Histogram vs PDF at tip of bar
figure
hiss = histogram(u_mc(:, end), 'Normalization', 'pdf');
set(hiss, 'DisplayName', 'Hist.')
hold on
plot_PDF_MC = plot(xmc_pdf, pdf_u_tip_mc, 'b-', 'LineWidth', 1);
set(plot_PDF_MC, 'DisplayName', 'PDF')
hl = legend('show', 'FontSize', 8, 'location', 'northeast');
set(hl, 'Interpreter', 'latex');
grid on
xlabel('$u(L)$', 'Interpreter', 'latex')
ylabel('$f\big(u(L)\big)$', 'Interpreter', 'latex')
xlim([15 27]); ylim([0 0.3]);
axis square
title('MC  histogram vs PDF at tip of bar')
hold off
%% Subfigure 3: CDF of MC vs PCE at tip of bar
figure
plot_CDF_MC = plot(xmc_cdf, cdf_u_tip_mc, 'b-', 'LineWidth', 1);
set(plot_CDF_MC, 'DisplayName', 'MC')
hold on
plot_CDF_PC = plot(xpc_cdf, cdf_u_tip_pc, 'r--', 'LineWidth', 1); % 'k-+' 'g-.' 'r--'
set(plot_CDF_PC, 'DisplayName', sprintf('PC'))
hl = legend('show', 'FontSize', 8, 'location', 'southeast');
set(hl, 'Interpreter', 'latex');
grid on
xlabel('$u(L)$', 'Interpreter', 'latex')
ylabel('$F\big(u(L)\big)$', 'Interpreter', 'latex')
axis square
title('CDF of MC vs PCE at tip of bar')
hold off
%% Subfigure 4: PDF of MC vs PCE at tip of bar
figure
plot_PDF_MC = plot(xmc_pdf, pdf_u_tip_mc, 'b-', 'LineWidth', 1);
set(plot_PDF_MC, 'DisplayName', 'MC')
hold on
plot_PDF_PC = plot(xpc_pdf, pdf_u_tip_pc, 'r--', 'LineWidth', 1); % 'k-+' 'g-.' 'r--'
set(plot_PDF_PC, 'DisplayName', sprintf('PC'))
hl = legend('show', 'FontSize', 8, 'location', 'northeast');
set(hl, 'Interpreter', 'latex');
grid on
xlabel('$u(L)$', 'Interpreter', 'latex')
ylabel('$f\big(u(L)\big)$', 'Interpreter', 'latex')
xlim([15 27]); ylim([0 0.3]);
axis square
title('PDF of MC vs PCE at tip of bar')
hold off
%% Heatmap of prior covaiance matrix
figure;
cm_viridis = viridis(100);
imagesc(C_u_pc)
sizeCov = size(C_u_pc, 1);
colormap(cm_viridis)
set(gca, 'xtick', [0:sizeCov / 2:sizeCov], 'ytick', [0:sizeCov / 2:sizeCov])
set(gca, 'XTickLabel', {'0', '50', '100'}, ...
  'YTickLabel', {'0', '50', '100'})
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$x$', 'Interpreter', 'latex')
axis square
colorbar
title('Prior Covariance Matrix')
%% Posterior dispalcement given experimental data.
figure;
title('Posterior displacement given experimental data')
plot_u_mean_pc = plot(nodeCoordinates, mu_u_pc, 'b-', 'LineWidth', 0.5);
set(plot_u_mean_pc, 'DisplayName', '$\mu_u$');
hold on
%
opts = {'EdgeColor', [0 0 1], 'FaceColor', [0 0 1], 'alpha', 0.5};
where = nodeCoordinates >- 1 & nodeCoordinates < nodeCoordinates(end) + 1;
plot_ci_c_u = fill_between(nodeCoordinates, (mu_u_pc - ci_u_pc), (mu_u_pc + ci_u_pc), where, opts{:});
set(plot_ci_c_u, 'DisplayName', '$95$\% CI')
hold on
% Create scatter plot for all points in one call
h = scatter(repmat(senCoor, 1, nrep), y_obs, 'filled', 'k', 'SizeData', 1, 'MarkerFaceAlpha', 0.5);
%
% Set legend only for one dataset and hide others
set(h, 'HandleVisibility', 'off'); % Hide all scatter plots from legend
set(h(1), 'DisplayName', 'obs. data', 'HandleVisibility', 'on'); % Show only first entry in legend
%
hold on
plot_u_y_mean = plot(nodeCoordinates, mu_u_y, 'r-', 'LineWidth', 0.5);
set(plot_u_y_mean, 'DisplayName', '$\mu_u|y$');
%
hold on
opts = {'EdgeColor', [1 0 0], 'FaceColor', [1 0 0], 'alpha', 0.3};
where = nodeCoordinates >- 1 & nodeCoordinates < nodeCoordinates(end) + 1;
plot_ci_u_y = fill_between(nodeCoordinates, (mu_u_y - ci_u_y), (mu_u_y + ci_u_y), where, opts{:});
set(plot_ci_u_y, 'DisplayName', '$95$\% CI$|y$')

%
hl = legend('show', 'FontSize', 8, 'location', 'northwest'); set(hl, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$u(x)$', 'Interpreter', 'latex');
ylim([-1 25]); xlim([-5 L + 5]);
axis square; grid on; hold off;
title('Posterior displacement given experimental data')
%% True system response given experimental data.
figure;
plot_u_mean_pc = plot(nodeCoordinates, mu_u_pc, 'b-', 'LineWidth', 0.5);
set(plot_u_mean_pc, 'DisplayName', '$\mu_u$');
hold on
%
opts = {'EdgeColor', [0 0 1], 'FaceColor', [0 0 1], 'alpha', 0.5};
where = nodeCoordinates >- 1 & nodeCoordinates < nodeCoordinates(end) + 1;
plot_ci_c_u = fill_between(nodeCoordinates, (mu_u_pc - ci_u_pc), (mu_u_pc + ci_u_pc), where, opts{:});
set(plot_ci_c_u, 'DisplayName', '$95$\% CI')
hold on
% Create scatter plot for all points in one call
h = scatter(repmat(senCoor, 1, nrep), y_obs, 'filled', 'k', 'SizeData', 1, 'MarkerFaceAlpha', 0.5);
%
% Set legend only for one dataset and hide others
set(h, 'HandleVisibility', 'off'); % Hide all scatter plots from legend
set(h(1), 'DisplayName', 'obs. data', 'HandleVisibility', 'on'); % Show only first entry in legend
%
plot_z = plot(senCoor, mu_z, 'k-', 'LineWidth', 1);
set(plot_z, 'DisplayName', '$z|y$')
hold on
%
opts = {'EdgeColor', [0 0 0], 'FaceColor', [0 0 0], 'alpha', 0.2};
where = senCoor >- 1 & senCoor < (senCoor(end) + 1);
plot_ci_z = fill_between(senCoor, (mu_z - ci_z), (mu_z + ci_z), where, opts{:});
set(plot_ci_z, 'DisplayName', '$95$\% CI of $z$')
%
hl = legend('show', 'FontSize', 8, 'location', 'northwest'); set(hl, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$u(x)$', 'Interpreter', 'latex');
ylim([-1 25]); xlim([-5 L + 5]);
axis square; grid on; hold off;
title('True system response given experimental data')
%% Heatmap of posterior displacement given experimental data
figure;
cm_viridis = viridis(100);
imagesc(C_u_y)
sizeCov = size(C_u_y, 1);
colormap(cm_viridis)
set(gca, 'xtick', [0:sizeCov / 2:sizeCov], 'ytick', [0:sizeCov / 2:sizeCov])
set(gca, 'XTickLabel', {'0', '50', '100'}, ...
  'YTickLabel', {'0', '50', '100'})
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$x$', 'Interpreter', 'latex')
axis square
colorbar
title('Posterior Covariance Matrix')
%% Heatmap of true system response given experimental data
figure;
cm_viridis = viridis(100);
imagesc(C_z)
sizeCov = size(C_z, 1);
colormap(cm_viridis)
set(gca, 'xtick', [0:sizeCov / 2:sizeCov], 'ytick', [0:sizeCov / 2:sizeCov])
set(gca, 'XTickLabel', {'0', '50', '100'}, ...
  'YTickLabel', {'0', '50', '100'})
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$x$', 'Interpreter', 'latex')
axis square
colorbar
title('True Covariance Matrix')
end
