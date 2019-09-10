function plot_statistic_log_v2(MESH,SETTINGS)
% Usage: plot_statistic_log_v2(MESH,SETTINGS)
%
% Purpose:
%   Plot nodes outside boundary, worst q, mean q, percetage under 0.6 and
%   mean misfit of the springs for each iteration
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Apr 2015
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if MESH.iter == 0; clf; end

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

subplot(3,1,1)
plot(MESH.iter,MESH.worst_q,'r+')
grid on
title('worst q')
hold on

subplot(3,1,2)
plot(MESH.iter,MESH.mean_q,'b+')
grid on
title('mean q')
hold on

subplot(3,1,3)
q_under_04        = MESH.q(MESH.q<0.4);
per_cent_under_04 = (size(q_under_04,1)/size(MESH.q,1))*100;
plot(MESH.iter, per_cent_under_04,'kx')
grid on
title('% q under 0.4')
hold on

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    saveas(gcf,[prefix 'Statistic_log_v2'],'fig');
end

end % END OF FUNCTION plot_statistic_log_v2