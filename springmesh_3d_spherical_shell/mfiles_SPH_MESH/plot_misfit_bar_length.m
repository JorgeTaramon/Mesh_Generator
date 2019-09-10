function plot_misfit_bar_length(MESH,SETTINGS)
% Usage: plot_misfit_bar_length(MESH,SETTINGS)
%
% Purpose:
%   Plot a histogram of how close each bar is to the desired bar length, L0
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Jan 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

rel_change_abs             = MESH.rel_change_abs;
mean_misfit_bar_length     = MESH.mean_misfit_bar_length;
mean_misfit_bar_length_tol = SETTINGS.mean_misfit_bar_length_tol;
x                          = (0.5:5:200);
nr                         = zeros(1,40);
N                          = size(rel_change_abs,1);
for i = 1:40
    r     = rel_change_abs(((i-1)/20<=rel_change_abs)&(rel_change_abs<i/20));
    nr(i) = 100*size(r,1)/N;
end
hold on
plot(x,nr,'b')
bar(mean_misfit_bar_length_tol*100,max(nr),'BarWidth',0.1,'EdgeColor',[1 0 0])
text(103,max(nr),...
    ['tolerance for mean misfit bar-length = ',sprintf('%.1f',100*mean_misfit_bar_length_tol),'%'],...
    'BackgroundColor',[1 1 1],'Margin',3,'EdgeColor','red');
bar(mean_misfit_bar_length*100,max(nr),'BarWidth',0.1)
text(130,max(nr)-2,...
    ['mean misfit bar-length = ',sprintf('%.1f',100*mean_misfit_bar_length),'%'],...
    'BackgroundColor',[1 1 1],'Margin',3,'EdgeColor','black');
title('Misfit bar-length')
xlabel('% of misfit')
ylabel('% of bars')
text(181,max(nr)-4,['iter = ',num2str(MESH.iter)],'BackgroundColor',[1 1 1],'Margin',2,'EdgeColor','black');

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Misfit_bar_length' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Misfit_bar_length' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_misfit_bar_length