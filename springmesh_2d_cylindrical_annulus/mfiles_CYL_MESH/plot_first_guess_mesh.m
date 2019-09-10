function plot_first_guess_mesh(FIRST_GUESS_MESH,SETTINGS)
% Usage: plot_first_guess_mesh(FIRST_GUESS_MESH,SETTINGS)
%
% Purpose:
%   Plot of first guess mesh
%
% Input:
%   FIRST_GUESS_MESH : [structure] : structure containing 1st guess for the
%                                    mesh
%   SETTINGS         : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

figure(2); clf

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

plot(FIRST_GUESS_MESH.pfix(:,1),FIRST_GUESS_MESH.pfix(:,2),'xr')
hold on
if SETTINGS.r_int > 0 && ~isempty(FIRST_GUESS_MESH.pbnd1)
    plot(FIRST_GUESS_MESH.pbnd1(:,1),FIRST_GUESS_MESH.pbnd1(:,2),'xb')
end
plot(FIRST_GUESS_MESH.pbnd2(:,1),FIRST_GUESS_MESH.pbnd2(:,2),'+b')
if ~isempty(FIRST_GUESS_MESH.pint)
    plot(FIRST_GUESS_MESH.pint(:,1),FIRST_GUESS_MESH.pint(:,2),'*g')
end
title('first guess for nodes')
xlabel('X (km)')
ylabel('Y (km)')
axis equal
if SETTINGS.r_int > 0 && ~isempty(FIRST_GUESS_MESH.pbnd1) && ~isempty(FIRST_GUESS_MESH.pint)
    legend('pfix','pbnd1','pbnd2','pint')
elseif SETTINGS.r_int > 0 && ~isempty(FIRST_GUESS_MESH.pbnd1) && isempty(FIRST_GUESS_MESH.pint)
    legend('pfix','pbnd1','pbnd2')
elseif SETTINGS.r_int > 0 && isempty(FIRST_GUESS_MESH.pbnd1) && ~isempty(FIRST_GUESS_MESH.pint)
    legend('pfix','pbnd2','pint')
elseif SETTINGS.r_int > 0 && isempty(FIRST_GUESS_MESH.pbnd1) && isempty(FIRST_GUESS_MESH.pint)
    legend('pfix','pbnd2')
end

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'First_guess_mesh']);
        case 'fig'
            saveas(gcf,[prefix 'First_guess_mesh'],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_first_guess_mesh