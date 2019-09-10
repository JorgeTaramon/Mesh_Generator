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
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

GCOORD = MESH.GCOORD';
EL2NOD = MESH.EL2NOD';
x      = reshape(GCOORD(1,EL2NOD),3,[]);
y      = reshape(GCOORD(2,EL2NOD),3,[]);
c      = repmat(MESH.q,1,3)';

patch(x,y,c,'FaceColor','flat')
axis equal
caxis([0.5 1])
colorbar('EastOutside')
xmin = SETTINGS.x0 - SETTINGS.length/2;
text(xmin*0.8,-200,['iter = ',num2str(MESH.iter)],'HorizontalAlignment', ...
    'center','BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
xlabel('X (km)')
ylabel('Y (km)')
title('Mesh quality factor')

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

end % END OF FUNCTION plot_mesh_quality