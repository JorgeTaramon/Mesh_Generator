function [MESH,SETTINGS] = first_guess_using_loaded_mesh_as_pfix(SETTINGS,GUIDE_MESH)
% Usage: [MESH,SETTINGS] = first_guess_using_loaded_mesh_as_pfix(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Load a first guess mesh as pfix and create the rest of nodes according
%   parameters in SETTINGS.
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings and
%                              the guide mesh itself in polar coordinates
%                              (theta,r)
%
% Output:
%   MESH       : [structure] : structure containing 1st guess for the mesh
%   SETTINGS   : [structure] : structure containing mesh settings
%
% JMT Oct 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

filename                = [pwd '/' SETTINGS.first_guess_file];
filename(filename=='\') = '/';
try
    tmp = load(filename);
catch
    error('Cannot find mesh file %s',filename);
end
if ~isfield(tmp,'MESH')
    error('Meshfile %s does not contain structure MESH.',filename);
end
MESH = tmp.MESH; clear tmp
pfix = [MESH.pfix; MESH.pbnd1; MESH.pbnd2; MESH.pint];

if max(max(pfix)) == SETTINGS.r_int
    pbnd1 = [];
    %======================================================================
    % CREATE BOUNDARY 2 NODES AND INTERIOR NODES
    %======================================================================
    if ~isempty(GUIDE_MESH) == 1
        [~,~,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH);
        pint        = interior_nodes(SETTINGS,GUIDE_MESH);
    else
        [~,~,pbnd2] = boundary_nodes_no_guide(SETTINGS);
        pint        = interior_nodes_no_guide(SETTINGS);
    end
    SETTINGS.r_int  = MESH.r_int;
elseif max(max(pfix)) == SETTINGS.r_ext
    error('need to be tested')
    pbnd2 = [];
    %======================================================================
    % CREATE BOUNDARY 1 NODES AND INTERIOR NODES
    %======================================================================
    if ~isempty(GUIDE_MESH) == 1
        [~,pbnd1,~] = boundary_nodes(SETTINGS,GUIDE_MESH);
        pint        = interior_nodes(SETTINGS,GUIDE_MESH);
    else
        [~,pbnd1,~] = boundary_nodes_no_guide(SETTINGS);
        pint        = interior_nodes_no_guide(SETTINGS);
    end
    SETTINGS.r_ext  = MESH.r_ext;
else
    %======================================================================
    % THE MESH LOADED CAN BE A SMALL BALL WITHIN THE CURRENT MESH
    %======================================================================
    if ~isempty(GUIDE_MESH) == 1
        [pfix_new,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH);
        [pint]                 = interior_nodes(SETTINGS,GUIDE_MESH);
        % Remove nodes (in a square) where the loaded mesh is placed
        xmax = max(pfix(:,1)) + GUIDE_MESH.l0_ref/2;
        xmin = min(pfix(:,1)) - GUIDE_MESH.l0_ref/2;
        ymax = max(pfix(:,2)) + GUIDE_MESH.l0_ref/2;
        ymin = min(pfix(:,2)) - GUIDE_MESH.l0_ref/2;
        pfix_new(pfix_new(:,1) > xmin & pfix_new(:,1) < xmax & ...
                 pfix_new(:,2) > ymin & pfix_new(:,2) < ymax,:) = [];
        if SETTINGS.r_int > 0
            pbnd1(pbnd1(:,1) > xmin & pbnd1(:,1) < xmax & ...
                  pbnd1(:,2) > ymin & pbnd1(:,2) < ymax,:) = [];
        end
        pbnd2(pbnd2(:,1) > xmin & pbnd2(:,1) < xmax & ...
              pbnd2(:,2) > ymin & pbnd2(:,2) < ymax,:) = [];
        pint(pint(:,1) > xmin & pint(:,1) < xmax & ...
             pint(:,2) > ymin & pint(:,2) < ymax,:) = [];
        pfix = [pfix_new; pfix];
        
%         % Remove nodes (in a circle) where the loaded mesh is placed
%         xc = 0;    % x-coordinate of the center of the blob
%         yc = 6221; %6121; % y-coordinate of the center of the blob
%         r  = max(pfix(:,1)) + GUIDE_MESH.l0_ref/2;  % radius of the blob
%         pfix_new(sqrt( (pfix_new(:,1)-xc).^2 + (pfix_new(:,2)-yc).^2 ) <= r,:) = [];
%         if SETTINGS.r_int > 0
%             pbnd1(sqrt( (pbnd1(:,1)-xc).^2 + (pbnd1(:,2)-yc).^2 ) <= r,:) = [];
%         end
%         pbnd2(sqrt( (pbnd2(:,1)-xc).^2 + (pbnd2(:,2)-yc).^2 ) <= r,:) = [];
%         pint(sqrt( (pint(:,1)-xc).^2 + (pint(:,2)-yc).^2 ) <= r,:) = [];
%         pfix = [pfix_new; pfix];
    else
        [pfix_new,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS);
        [pint]                 = interior_nodes_no_guide(SETTINGS);
%         % Remove nodes (in a square) where the loaded mesh is placed
%         xmax = max(pfix(:,1)) + SETTINGS.h0/2;
%         xmin = min(pfix(:,1)) - SETTINGS.h0/2;
%         ymax = max(pfix(:,2)) + SETTINGS.h0/2;
%         ymin = min(pfix(:,2)) - SETTINGS.h0/2;
%         pfix_new(pfix_new(:,1) > xmin & pfix_new(:,1) < xmax & ...
%                  pfix_new(:,2) > ymin & pfix_new(:,2) < ymax,:) = [];
%         if SETTINGS.r_int > 0
%             pbnd1(pbnd1(:,1) > xmin & pbnd1(:,1) < xmax & ...
%                   pbnd1(:,2) > ymin & pbnd1(:,2) < ymax,:) = [];
%         end
%         pbnd2(pbnd2(:,1) > xmin & pbnd2(:,1) < xmax & ...
%               pbnd2(:,2) > ymin & pbnd2(:,2) < ymax,:) = [];
%         pint(pint(:,1) > xmin & pint(:,1) < xmax & ...
%              pint(:,2) > ymin & pint(:,2) < ymax,:) = [];
%         pfix = [pfix_new; pfix];
        
        % Remove nodes (in a circle) where the loaded mesh is placed
        xc = 0;    % x-coordinate of the center of the blob
        yc = 6271; % y-coordinate of the center of the blob
        r  = 27;  % radius of the blob
        pfix_new(sqrt( (pfix_new(:,1)-xc).^2 + (pfix_new(:,2)-yc).^2 ) <= r,:) = [];
        if SETTINGS.r_int > 0
            pbnd1(sqrt( (pbnd1(:,1)-xc).^2 + (pbnd1(:,2)-yc).^2 ) <= r,:) = [];
        end
        pbnd2(sqrt( (pbnd2(:,1)-xc).^2 + (pbnd2(:,2)-yc).^2 ) <= r,:) = [];
        pint(sqrt( (pint(:,1)-xc).^2 + (pint(:,2)-yc).^2 ) <= r,:) = [];
        pfix = [pfix_new; pfix];
    end
end

%==========================================================================
% OUTPUT DATA
%==========================================================================
GCOORD     = [pfix;pbnd1;pbnd2;pint];
EL2NOD     = delaunay(GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_POL = cartesian2polar(GCOORD);
EL2NOD     = EL2NOD(~(sum(ismember(EL2NOD,find(abs(GCOORD_POL(:,2)-SETTINGS.r_int) < 1e-8)),2)==3),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_POL(:,2)-SETTINGS.r_int) < 1e-8);
% elements_with_3_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 3;
% remove_elements_with_3_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_3_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------
MESH.pfix   = pfix;
MESH.pbnd1  = pbnd1;
MESH.pbnd2  = pbnd2;
MESH.pint   = pint;
MESH.GCOORD = GCOORD;
MESH.EL2NOD = EL2NOD;
MESH.q      = mesh_quality(MESH.GCOORD,MESH.EL2NOD);

%==========================================================================
% ADD pfix refined region
%==========================================================================
if strcmp(SETTINGS.refinement,'guide_mesh') && strcmp(SETTINGS.bnds_at_ref_region,'fixed')
    MESH = pfix_refined_region(MESH,GUIDE_MESH,SETTINGS);
end

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    plot_first_guess_mesh(MESH,SETTINGS)
    MESH.iter      = 0;
    SETTINGS.iplot = 0;
    figure(15)
    plot_mesh(MESH,SETTINGS)
    figure(16)
    MESH.worst_q = min(MESH.q);
    MESH.mean_q  = mean(MESH.q);
    plot_statistic_log_v2(MESH,SETTINGS)
end

end % END OF FUNCTION first_guess_using_loaded_mesh_as_pfix