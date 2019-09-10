function plot_guide_mesh(GUIDE_MESH,SETTINGS)
% Usage: plot_guide_mesh(GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Plot of guide mesh
%
% Input:
%   GUIDE_MESH : [structure] : structure containing guide mesh
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   plot
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

figure(1); clf

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

plot(GUIDE_MESH.p_coarse_guide(:,1),GUIDE_MESH.p_coarse_guide(:,2),'ob')
hold on
plot(GUIDE_MESH.p_tran_guide(:,1),GUIDE_MESH.p_tran_guide(:,2),'og')
plot(GUIDE_MESH.p_ref_guide(:,1),GUIDE_MESH.p_ref_guide(:,2),'or')
axis([-25000 25000 -3500 500])
legend(strcat('Coarse Zone (l_0 = ',num2str(GUIDE_MESH.l0_coarse),' km)')               , ...
       strcat('Transition Zone (boundary) (l_0 = ',num2str(GUIDE_MESH.l0_coarse),' km)'), ...
       strcat('Refined Zone (l_0 = ',num2str(GUIDE_MESH.l0_ref),' km)')                 , ...
       'Location','SouthOutside')
trimesh(GUIDE_MESH.EL2NOD_GUIDE,GUIDE_MESH.GCOORD_GUIDE(:,1),GUIDE_MESH.GCOORD_GUIDE(:,2),'Color',[0 0 0]);
xlabel('X (km)')
ylabel('Z (km)')
title('Guide Mesh')
text(GUIDE_MESH.GCOORD_GUIDE(:,1)+3,...
     GUIDE_MESH.GCOORD_GUIDE(:,2)+5,...
     num2str(GUIDE_MESH.L0_guide,4),'Color',[0 0 0],'Fontsize',10)

if SETTINGS.save_figs
    prefix  = [SETTINGS.outdir '/'];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Guide_mesh']);
        case 'fig'
            saveas(gcf,[prefix 'Guide_mesh'],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
end

end % END OF FUNCTION plot_guide_mesh