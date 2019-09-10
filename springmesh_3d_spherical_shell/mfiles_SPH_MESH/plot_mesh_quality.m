function plot_mesh_quality(MESH,SETTINGS)
% Usage: plot_mesh_quality(MESH,SETTINGS)
%
% Purpose:
%   Plot mesh quality
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Jul 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

figure(12);clf
if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end
GCOORD = MESH.GCOORD';
EL2NOD = MESH.EL2NOD';
q      = MESH.q;
x      = reshape(GCOORD(1,EL2NOD),4,[]);
y      = reshape(GCOORD(2,EL2NOD),4,[]);
z      = reshape(GCOORD(3,EL2NOD),4,[]);
c      = repmat(q,1,4)';

patch(x,y,z,c,'FaceColor','flat','FaceAlpha',0.4)
caxis([SETTINGS.q_tol 1])
colorbar('EastOutside')
grid on
view(142.5,30)
axis equal
xlabel('X');ylabel('Y');zlabel('Z');
title('Mesh quality')

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Mesh_quality' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Mesh_quality' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

figure(13);clf
if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end
simpplot(MESH.GCOORD,MESH.EL2NOD,'p(:,1)<0')
view(142.5,30)
xlabel('X');ylabel('Y');zlabel('Z');
title('Mesh cross-section')
if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Mesh_cross_section' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Mesh_cross_section' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_mesh_quality