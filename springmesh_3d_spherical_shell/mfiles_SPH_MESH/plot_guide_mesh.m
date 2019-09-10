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
% JMT Jun 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

figure(1); clf

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

scatter3(GUIDE_MESH.p_coarse_guide(:,1),...
         GUIDE_MESH.p_coarse_guide(:,2),...
         GUIDE_MESH.p_coarse_guide(:,3)/10,100,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
hold on
scatter3(GUIDE_MESH.p_tran_guide(:,1),...
         GUIDE_MESH.p_tran_guide(:,2),...
         GUIDE_MESH.p_tran_guide(:,3)/10,100,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
scatter3(GUIDE_MESH.p_ref_guide(:,1),...
         GUIDE_MESH.p_ref_guide(:,2),...
         GUIDE_MESH.p_ref_guide(:,3)/10,100,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
axis([0 180 0 360 347.1 637.1])
set(gca,'XTick',[0 30 60 90 120 150 180],...
        'YTick',[0 30 60 90 120 150 180 210 240 270 300 330 360],...
        'ZTick',[347.1 400.0 450.0 500.0 550.0 600.0 637.1])
view(120,45)
legend(strcat('Coarse Zone (l_0 = ',num2str(GUIDE_MESH.l0_coarse),' km)'),...
       strcat('Transition Zone (boundary) (l_0 = ',num2str(GUIDE_MESH.l0_coarse),' km)'),...
       strcat('Refined Zone (l_0 = ',num2str(GUIDE_MESH.l0_ref),' km)'),...
       'Location','SouthOutside')
GCOORD_GUIDE_scaled      = GUIDE_MESH.GCOORD_GUIDE;
GCOORD_GUIDE_scaled(:,3) = GCOORD_GUIDE_scaled(:,3)/10;
tetramesh(GUIDE_MESH.EL2NOD_GUIDE,GCOORD_GUIDE_scaled,'FaceColor',[0 0 0],'FaceAlpha',0);  
xlabel('\theta ({\circ})');
ylabel('\phi ({\circ})');
zlabel('r / 10 (km)')
title('Guide Mesh')
% text(GUIDE_MESH.GCOORD_GUIDE(:,1)+3,...
%      GUIDE_MESH.GCOORD_GUIDE(:,2)+3,...
%      GUIDE_MESH.GCOORD_GUIDE(:,3)/10+5,...
%      num2str(GUIDE_MESH.L0_guide,4),'Color',[0 0 0],'Fontsize',10)

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