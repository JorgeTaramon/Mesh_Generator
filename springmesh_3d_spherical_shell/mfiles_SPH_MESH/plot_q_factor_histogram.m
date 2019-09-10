function plot_q_factor_histogram(q,iter,SETTINGS)
% Usage: plot_q_factor_histogram(q,iter,SETTINGS)
%
% Purpose:
%   Plot a histogram for the quality factor
%
% Input:
%   q        : [column vector] : quality factor for each tetrahedron
%   iter     : [scalar]        : iteration number
%   SETTINGS : [structure]     : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Apr 2015
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

q1   = q((0<q)&(q<=0.1));   % 0 < q <= 0.1
nq1  = size(q1,1);
q2   = q((0.1<q)&(q<=0.2)); % 0.1 < q <= 0.2 
nq2  = size(q2,1);
q3   = q((0.2<q)&(q<=0.3)); % 0.2 < q <= 0.3
nq3  = size(q3,1);
q4   = q((0.3<q)&(q<=0.4)); % 0.3 < q <= 0.4
nq4  = size(q4,1);
q5   = q((0.4<q)&(q<=0.5)); % 0.4 < q <= 0.5
nq5  = size(q5,1);
q6   = q((0.5<q)&(q<=0.6)); % 0.5 < q <= 0.6
nq6  = size(q6,1);
q7   = q((0.6<q)&(q<=0.7)); % 0.6 < q <= 0.7
nq7  = size(q7,1);
q8   = q((0.7<q)&(q<=0.8)); % 0.7 < q <= 0.8
nq8  = size(q8,1);
q9   = q((0.8<q)&(q<=0.9)); % 0.8 < q <= 0.9
nq9  = size(q9,1);
q10  = q((0.9<q)&(q<=1));   % 0.9 < q <= 1
nq10 = size(q10,1);
x    =[0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95];
nq   = [nq1 nq2 nq3 nq4 nq5 nq6 nq7 nq8 nq9 nq10];
bar(x,nq,'BarWidth',1,'FaceColor',[.7 .9 .7])
set(gca, 'YScale', 'log')
set(gca,'YGrid','on')
title('Histogram of quality factor, q')
xlabel('quality factor')
text(0.1, 200,['worst q = ',sprintf('%.3f',min(q))],...
    'BackgroundColor',[1 1 1],'Margin',5,'LineWidth',2,'EdgeColor','red');
text(0.1,500,['iter = ',num2str(iter)],...
    'BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Hist_q_factor' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Hist_q_factor' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

% =========================================================================
% Alternative and easier way to plot the histogram:
%
% hist(q,[0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]) 
%
% The problem is that the log scale doesn't work as in the previous figure.
% =========================================================================

end % END OF FUNCTION plot_q_factor_histogram