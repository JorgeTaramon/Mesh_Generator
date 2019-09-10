function plot_net_change_nodes(MESH,SETTINGS)
% Usage: plot_net_change_nodes(MESH,SETTINGS)
%
% Purpose:
%   Plot fraction of nodes changed and net change of nodes
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

if MESH.add_reject_counter == 1; clf; end

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

if SETTINGS.r_int > 0 && ~isempty(MESH.pbnd1)
    subplot(3,1,1)
    bar(MESH.add_reject_counter,MESH.fraction_bnd1_changed,'BarWidth',1,'FaceColor',[0 0 1]);
    hold on
    bar(MESH.add_reject_counter,MESH.net_change_bnd1,'BarWidth',1,'FaceColor',[.7 .9 .7]);
    ylim([0 1])
    legend('Fraction of changed bnd1 nodes','Net change of bnd1 nodes')
end
subplot(3,1,2)
bar(MESH.add_reject_counter,MESH.fraction_bnd2_changed,'BarWidth',1,'FaceColor',[0 0 1]);
hold on
bar(MESH.add_reject_counter,MESH.net_change_bnd2,'BarWidth',1,'FaceColor',[.7 .9 .7]);
ylim([0 1])
legend('Fraction of changed bnd2 nodes','Net change of bnd2 nodes')
subplot(3,1,3)
bar(MESH.add_reject_counter,MESH.fraction_int_changed,'BarWidth',1,'FaceColor',[0 0 1]);
hold on
bar(MESH.add_reject_counter,MESH.net_change_int,'BarWidth',1,'FaceColor',[.7 .9 .7]);
ylim([0 1])
legend('Fraction of changed int nodes','Net change of int nodes')

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Net_nodes_change']);
        case 'fig'
            saveas(gcf,[prefix 'Net_nodes_change'],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_net_change_nodes