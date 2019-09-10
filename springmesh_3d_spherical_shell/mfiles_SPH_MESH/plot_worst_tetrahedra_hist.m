function plot_worst_tetrahedra_hist(q,iter,SETTINGS)
% Usage: plot_worst_tetrahedra_hist(q,iter,SETTINGS)
%
% Purpose:
%   Plot a histogram for the quality factor of the worst tetrahedra
%
% Input:
%   q        : [column vector] : quality factor for each tetrahedron
%   iter     : [scalar]        : iteration number
%   SETTINGS : [structure]     : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Jul 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

ibad{1}  = find(q<1e-1);
ibad{2}  = find(q>1e-1 & q<2e-1);
ibad{3}  = find(q>2e-1 & q<3e-1);
ibad{4}  = find(q>3e-1 & q<4e-1);

q0_q1 = length(ibad{1});
q1_q2 = length(ibad{2});
q2_q3 = length(ibad{3});
q3_q4 = length(ibad{4});
x     = [0.05, 0.15, 0.25, 0.35];
nr    = [q0_q1 q1_q2 q2_q3 q3_q4];
bar(x,nr,'BarWidth',1,'FaceColor',[.7 .9 .7])
set(gca,'YGrid','on')
text(0.015,q0_q1+1,[      'q < 0.1 (',num2str(q0_q1),'/',num2str(length(q)),')'],'BackgroundColor',[0.6875 0.8750 0.8984],'Margin',3,'EdgeColor','black');
text(0.1  ,q1_q2+1,['0.1 < q < 0.2 (',num2str(q1_q2),'/',num2str(length(q)),')'],'BackgroundColor',[0.6875 0.8750 0.8984],'Margin',3,'EdgeColor','black');
text(0.2  ,q2_q3+1,['0.2 < q < 0.3 (',num2str(q2_q3),'/',num2str(length(q)),')'],'BackgroundColor',[0.6875 0.8750 0.8984],'Margin',3,'EdgeColor','black');
text(0.3  ,q3_q4-2,['0.3 < q < 0.4 (',num2str(q3_q4),'/',num2str(length(q)),')'],'BackgroundColor',[0.6875 0.8750 0.8984],'Margin',3,'EdgeColor','black');
text(0.015,q3_q4-7,['worst q = ',sprintf('%.3f',min(q))],'BackgroundColor',[1 1 1],'Margin',5,'LineWidth',2,'EdgeColor','red');
text(0.015,q3_q4-2,['iter = ',num2str(iter)],'BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black');

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Hist_worst_q_factor' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Hist_worst_q_factor' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_worst_tetrahedra_hist