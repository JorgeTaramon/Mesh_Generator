function MESH = plot_output_data(MESH,SETTINGS)
% Usage: MESH = plot_output_data(MESH,SETTINGS)
%
% Purpose:
%   Plot different output data
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   MESH     : [structure] : structure containing the mesh
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

% Plot nodes outside boundary, worst q, mean q, % of q under 0.6 and mean misfit of springs
figure(6)
plot_statistic_log(MESH,SETTINGS)

% Plot a histogram for the quality factor
figure(7); clf
plot_percentage_q_factor(MESH.q,MESH.iter,SETTINGS,MESH.s)
% plot_q_factor_histogram(MESH.q,MESH.iter,SETTINGS)

% Plot cumulative percentage of quality factor
figure(8); clf
plot_cumulative_percentage_q_factor(MESH.q,MESH.s,MESH.iter,SETTINGS)

% Plot a histogram of percentage of misfit for the bars
figure(9); clf
plot_misfit_bar_length(MESH,SETTINGS)

% Plot hisogram for the worst tetrahedra (q < 0.4)
figure(19); clf
plot_worst_tetrahedra_hist(MESH.q,MESH.iter,SETTINGS)

% Plot 1st iteration where worst_q reaches 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55 and 0.60
figure(11); clf
MESH = plot_first_iter_reaching_different_worst_q_values(MESH,SETTINGS);

% Plot the mesh
figure(15); clf
plot_mesh(MESH,SETTINGS)

% Plot worst q, mean q and % of q under 0.6
figure(16)
plot_statistic_log_v2(MESH,SETTINGS)

end % END OF FUNCTION plot_output_data