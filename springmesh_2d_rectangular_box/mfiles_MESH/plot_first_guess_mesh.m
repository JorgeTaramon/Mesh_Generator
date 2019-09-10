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
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

figure(2); clf

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

plot(FIRST_GUESS_MESH.pfix(:,1),FIRST_GUESS_MESH.pfix(:,2),'xr')
hold on
plot(FIRST_GUESS_MESH.pbnd(:,1),FIRST_GUESS_MESH.pbnd(:,2),'xb')
plot(FIRST_GUESS_MESH.pint(:,1),FIRST_GUESS_MESH.pint(:,2),'*g')
title('first guess for nodes')
xlabel('X (km)')
ylabel('Y (km)')
% axis equal
legend('pfix','pbnd','pint')

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