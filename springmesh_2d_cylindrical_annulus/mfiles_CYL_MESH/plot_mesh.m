function plot_mesh(MESH,SETTINGS)
% Usage: plot_mesh(MESH,SETTINGS)
%
% Purpose:
%   Plot mesh
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Jun 2017
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
c      = repmat(0.90*[1 1 1],size(MESH.q,1),1)';
patch(x,y,c,'FaceColor',0.90*[1 1 1])
axis equal
grid on
axis([-SETTINGS.r_ext SETTINGS.r_ext -SETTINGS.r_ext SETTINGS.r_ext])
text(-SETTINGS.r_ext*0.9,SETTINGS.r_ext*0.9,['iter = ',num2str(MESH.iter)],...
    'BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
xlabel('X (km)')
ylabel('Y (km)')
title('Mesh')

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    saveas(gcf,[prefix 'Mesh' suffix],'fig');
end

end % END OF FUNCTION plot_mesh