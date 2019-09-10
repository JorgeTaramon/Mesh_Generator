function [MESH,SETTINGS] = first_guess_using_loaded_mesh_as_pfix(SETTINGS,GUIDE_MESH)
% Usage: [MESH,SETTINGS] = first_guess_using_loaded_mesh_as_pfix(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Load a first guess mesh as pfix and create the rest of nodes according
%   with the radius chosen in SETTINGS. This is useful to create spherical
%   shells with interior layers
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings and
%                              the guide mesh itself in spherical
%                              coordinates (theta,phi,r)
%
% Output:
%   MESH       : [structure] : structure containing 1st guess for the mesh
%   SETTINGS   : [structure] : structure containing mesh settings
%
% JMT Jan 2017
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
    %==========================================================================
    % CREATE BOUNDARY 2 NODES AND INTERIOR NODES
    %==========================================================================
    if ~isempty(GUIDE_MESH) == 1
        [~,~,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH);
        pint        = interior_nodes(SETTINGS,GUIDE_MESH);
    else
        [~,~,pbnd2] = boundary_nodes_no_guide(SETTINGS);
        pint        = interior_nodes_no_guide(SETTINGS);
    end
    SETTINGS.r_int  = MESH.r_int;
elseif max(max(pfix)) == SETTINGS.r_ext
%     error('need to be tested')
%     pbnd2 = [];
%     %==========================================================================
%     % CREATE BOUNDARY 1 NODES AND INTERIOR NODES
%     %==========================================================================
%     if ~isempty(GUIDE_MESH) == 1
%         [~,pbnd1,~] = boundary_nodes(SETTINGS,GUIDE_MESH);
%         pint        = interior_nodes(SETTINGS,GUIDE_MESH);
%     else
%         [~,pbnd1,~] = boundary_nodes_no_guide(SETTINGS);
%         pint        = interior_nodes_no_guide(SETTINGS);
%     end
%     SETTINGS.r_ext  = MESH.r_ext;
else
    %======================================================================
    % THE MESH LOADED CAN BE A SMALL BALL WITHIN THE CURRENT MESH
    %======================================================================
    if ~isempty(GUIDE_MESH) == 1
        [pfix_new,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH);
        [pint]                 = interior_nodes(SETTINGS,GUIDE_MESH);
    else
        [pfix_new,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS);
        [pint]                 = interior_nodes_no_guide(SETTINGS);
    end
    % Remove nodes (in a sphere) where the loaded mesh is placed
    xc = 0;    % x-coordinate of the center of the sphere
    yc = 6271; % y-coordinate of the center of the sphere
    zc = 0;    % z-coordinate of the center of the sphere
    r  = 85; %GUIDE_MESH.l0_ref/2;  % radius of the sphere
    pfix_new(sqrt( (pfix_new(:,1)-xc).^2 + ...
                   (pfix_new(:,2)-yc).^2 + ...
                   (pfix_new(:,3)-zc).^2 ) <= r,:) = [];
    if SETTINGS.r_int > 0
        pbnd1(sqrt( (pbnd1(:,1)-xc).^2 + ...
                    (pbnd1(:,2)-yc).^2 + ...
                    (pbnd1(:,3)-zc).^2) <= r,:) = [];
    end
    pbnd2(sqrt( (pbnd2(:,1)-xc).^2 + ...
                (pbnd2(:,2)-yc).^2 + ...
                (pbnd2(:,3)-zc).^2 ) <= r,:) = [];
    pint(sqrt( (pint(:,1)-xc).^2 + ...
               (pint(:,2)-yc).^2 + ...
               (pint(:,3)-zc).^2 ) <= r,:) = [];
    pfix = [pfix_new; pfix];
end

%==========================================================================
% OUTPUT DATA
%==========================================================================
GCOORD     = [pfix;pbnd1;pbnd2;pint];
EL2NOD     = delaunay(GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH = cartesian2spherical(GCOORD);
EL2NOD     = EL2NOD(~(sum(ismember(EL2NOD,find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8);
% elements_with_4_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 4;
% remove_elements_with_4_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_4_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------
MESH.pfix   = pfix;
MESH.pbnd1  = pbnd1;
MESH.pbnd2  = pbnd2;
MESH.pint   = pint;
MESH.GCOORD = GCOORD;
MESH.EL2NOD = EL2NOD;
MESH.q      = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
MESH.s      = shape_measure(MESH.GCOORD,MESH.EL2NOD);
MESH.r_int  = SETTINGS.r_int;
MESH.r_ext  = SETTINGS.r_ext;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    plot_first_guess_mesh(MESH,SETTINGS)
end

end % END OF FUNCTION first_guess_using_loaded_mesh_as_pfix