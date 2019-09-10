function plot_statistic_log(MESH,SETTINGS)
% Usage: plot_statistic_log(MESH,SETTINGS)
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

if MESH.iter == 1; clf; end

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

subplot(5,1,1)
plot(MESH.iter,MESH.n_nod_out,'g+')
title('nodes outside boundary')
hold on

subplot(5,1,2)
plot(MESH.iter,MESH.worst_q,'r+')
title('worst q')
hold on

subplot(5,1,3)
plot(MESH.iter,MESH.mean_q,'b+')
title('mean q')
hold on

subplot(5,1,4)
q_under_06        = MESH.q(MESH.q<0.6);
per_cent_under_06 = (size(q_under_06,1)/size(MESH.q,1))*100;
plot(MESH.iter, per_cent_under_06,'kx')
title('% q under 0.6')
hold on

subplot(5,1,5)
plot(MESH.iter,MESH.mean_misfit_bar_length,'k+')
title('mean misfit of the springs')
hold on

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Statistic_log']);
        case 'fig'
            saveas(gcf,[prefix 'Statistic_log'],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_statistic_log