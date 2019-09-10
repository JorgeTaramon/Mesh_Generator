function plot_output_data(MESH,SETTINGS)
% Usage: plot_output_data(MESH,SETTINGS)
%
% Purpose:
%   Plot different output data
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

% Plot nodes outside boundary, worst q, mean q, % of q under 0.6 and mean misfit of springs
figure(10)
plot_statistic_log(MESH,SETTINGS)

% Plot a histogram for the quality factor
figure(11); clf
% plot_percentage_q_factor(MESH.q,MESH.iter,SETTINGS)
plot_q_factor_histogram(MESH.q,MESH.iter,SETTINGS)

% Plot cumulative percentage of quality factor
figure(12); clf
plot_cumulative_percentage_q_factor(MESH.q,MESH.iter,SETTINGS)

% Plot a histogram of percentage of misfit for the bars
figure(13); clf
plot_misfit_bar_length(MESH,SETTINGS)

% Plot the mesh quality
figure(14); clf
plot_mesh_quality(MESH,SETTINGS)

% Plot the mesh
figure(15); clf
plot_mesh(MESH,SETTINGS)

% Plot worst q, mean q and % of q under 0.6
figure(16)
plot_statistic_log_v2(MESH,SETTINGS)

end % END OF FUNCTION plot_output_data