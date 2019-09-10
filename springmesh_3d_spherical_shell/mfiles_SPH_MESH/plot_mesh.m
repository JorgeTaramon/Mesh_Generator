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
r_ext = SETTINGS.r_ext;
simpplot(MESH.GCOORD,MESH.EL2NOD,'p(:,1)<0')
view(142.5,30)
xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)')
axis equal
axis([-r_ext r_ext -r_ext r_ext -r_ext r_ext])
grid on
text(r_ext/2,-r_ext,1.5*r_ext,['iter = ',num2str(MESH.iter)], ...
    'BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(MESH.iter,5)];
    saveas(gcf,[prefix 'Mesh_cross_section' suffix],'fig');
end

end % END OF FUNCTION plot_mesh


