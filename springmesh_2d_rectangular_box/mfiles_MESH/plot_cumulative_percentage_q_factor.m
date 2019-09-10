function plot_cumulative_percentage_q_factor(q,iter,SETTINGS)
% Usage: plot_cumulative_percentage_q_factor(q,iter,SETTINGS)
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

[nq,x] = compute_pct(q);

plot(x,nq,'r')
hold on
axis([0 1 0 100])
grid on
title('Cumulative percentage for quality factor')
xlabel('quality factor')
ylabel('cumulative % of quality factor')
text(0.035, 86,['worst q = ',sprintf('%.3f',min(q))],...
    'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor','black');
text(0.035,79,['iter = ',num2str(iter)],...
    'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor','black');

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Cumulative_percentage_q_factor' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Cumulative_percentage_q_factor' suffix],'fig');
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

function [n,x] = compute_pct(q)

N    = size(q,1);
q1   = q((0<q)&(q<=0.05));    
nq1  = 100*size(q1,1)/N;         % percentage of tetrahedra with 0 < q <= 0.05
q2   = q((0.05<q)&(q<=0.10));
nq2  = nq1  + 100*size(q2,1)/N;  % percentage of tetrahedra with 0 < q <= 0.10 
q3   = q((0.10<q)&(q<=0.15));
nq3  = nq2  + 100*size(q3,1)/N;  % percentage of tetrahedra with 0 < q <= 0.15
q4   = q((0.15<q)&(q<=0.20));
nq4  = nq3  + 100*size(q4,1)/N;  % percentage of tetrahedra with 0 < q <= 0.20
q5   = q((0.20<q)&(q<=0.25));
nq5  = nq4  + 100*size(q5,1)/N;  % percentage of tetrahedra with 0 < q <= 0.25
q6   = q((0.25<q)&(q<=0.30));
nq6  = nq5  + 100*size(q6,1)/N;  % percentage of tetrahedra with 0 < q <= 0.30
q7   = q((0.30<q)&(q<=0.35));
nq7  = nq6  + 100*size(q7,1)/N;  % percentage of tetrahedra with 0 < q <= 0.35
q8   = q((0.35<q)&(q<=0.40));
nq8  = nq7  + 100*size(q8,1)/N;  % percentage of tetrahedra with 0 < q <= 0.40
q9   = q((0.40<q)&(q<=0.45));
nq9  = nq8  + 100*size(q9,1)/N;  % percentage of tetrahedra with 0 < q <= 0.45
q10  = q((0.45<q)&(q<=0.50));
nq10 = nq9  +100*size(q10,1)/N;  % percentage of tetrahedra with 0 < q <= 0.50
q11  = q((0.50<q)&(q<=0.55));
nq11 = nq10 + 100*size(q11,1)/N; % percentage of tetrahedra with 0 < q <= 0.55
q12  = q((0.55<q)&(q<=0.60));
nq12 = nq11 + 100*size(q12,1)/N; % percentage of tetrahedra with 0 < q <= 0.60 
q13  = q((0.60<q)&(q<=0.65));
nq13 = nq12 + 100*size(q13,1)/N; % percentage of tetrahedra with 0 < q <= 0.65
q14  = q((0.65<q)&(q<=0.70));
nq14 = nq13 + 100*size(q14,1)/N; % percentage of tetrahedra with 0 < q <= 0.70
q15  = q((0.70<q)&(q<=0.75));
nq15 = nq14 + 100*size(q15,1)/N; % percentage of tetrahedra with 0 < q <= 0.75
q16  = q((0.75<q)&(q<=0.80));
nq16 = nq15 + 100*size(q16,1)/N; % percentage of tetrahedra with 0 < q <= 0.80
q17  = q((0.80<q)&(q<=0.85));
nq17 = nq16 + 100*size(q17,1)/N; % percentage of tetrahedra with 0 < q <= 0.85
q18  = q((0.85<q)&(q<=0.90));
nq18 = nq17 + 100*size(q18,1)/N; % percentage of tetrahedra with 0 < q <= 0.90
q19  = q((0.90<q)&(q<=0.95));
nq19 = nq18 + 100*size(q19,1)/N; % percentage of tetrahedra with 0 < q <= 0.95
q20  = q((0.95<q)&(q<=1));
nq20 = nq19 + 100*size(q20,1)/N; % percentage of tetrahedra with 0 < q <= 1
x    = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ...
        0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1];
n    = [nq1  nq2  nq3  nq4  nq5  nq6  nq7  nq8  nq9  nq10 ...
        nq11 nq12 nq13 nq14 nq15 nq16 nq17 nq18 nq19 nq20];
end