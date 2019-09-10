function MESH = pfix_refined_region(MESH,GUIDE_MESH,SETTINGS)
% Usage: MESH = pfix_refined_region(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Generation of pfix for embedded high resolution sub-region
%
% Input:
%   MESH       : [structure] : structure containing the mesh
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   MESH       : [structure] : structure containing the mesh
%
% JMT Jul 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
deg2rad      = pi/180;
r_ext        = SETTINGS.r_ext;      % outer radius (km) of cylindrical annulus
theta0       = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
d_ref        = GUIDE_MESH.d_ref;    % refined zone depth (km)
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees
theta_ref_l  = theta0-w_ref_deg/2;  % colatitude of the left boundary in the refined zone
theta_ref_r  = theta0+w_ref_deg/2;  % colatitude of the right boundary in the refined zone
l0_ref       = GUIDE_MESH.l0_ref;   % desired spring length (km) for refined zone
pfix         = MESH.pfix;
pbnd2        = MESH.pbnd2;

%===============================================================================
% REMOVE PFIX AND PBND ON BOUNDARY 2 (OUTER BOUNDARY) INSIDE THE REFINED REGION
%===============================================================================
pfix_pol     = cartesian2polar(pfix);
pfix(pfix_pol(:,1) >= theta_ref_l & pfix_pol(:,1) <= theta_ref_r & pfix_pol(:,2) > r_ext - 1e-8,:) = [];

pbnd2_pol    = cartesian2polar(pbnd2);
pbnd2(pbnd2_pol(:,1) >= theta_ref_l & pbnd2_pol(:,1) <= theta_ref_r,:) = [];

%==========================================================================
% DEFINE PFIX REFINED REGION
%==========================================================================
% pfix vertices refined region
pfix_vertex_pol    = [theta_ref_l  r_ext - d_ref ; ...
                      theta_ref_r  r_ext - d_ref ; ...
                      theta_ref_r  r_ext         ; ...
                      theta_ref_l  r_ext ]; % vertex points COUNTER-CLOCKWISE
pfix_vertex        = polar2cartesian(pfix_vertex_pol);

% pfix bottom ref region
l0_ref_degree_bot  = l0_ref/((r_ext - d_ref)*deg2rad);
arc_length_bot     = abs((theta_ref_l + l0_ref_degree_bot) - (theta_ref_r - l0_ref_degree_bot));
theta_bot          = linspace(theta_ref_l + l0_ref_degree_bot, theta_ref_r - l0_ref_degree_bot, round(arc_length_bot/l0_ref_degree_bot))';
r_bot              = (r_ext - d_ref)*ones(size(theta_bot,1),1);
p_fix_ref_bot_pol  = [theta_bot r_bot];
p_fix_ref_bot      = polar2cartesian(p_fix_ref_bot_pol); 

% pfix right side ref region
r_length_right     = abs((r_ext - d_ref + l0_ref) - (r_ext - l0_ref));
r_right            = linspace(r_ext - d_ref + l0_ref, r_ext - l0_ref, round(r_length_right/l0_ref))';
theta_right        = theta_ref_r*ones(size(r_right,1),1);
pfix_ref_right_pol = [theta_right r_right];
pfix_ref_right     = polar2cartesian(pfix_ref_right_pol);

% pfix top ref region
l0_ref_degree_top  = l0_ref/(r_ext*deg2rad);
arc_length_top     = abs((theta_ref_r - l0_ref_degree_top) - (theta_ref_l + l0_ref_degree_top));
theta_top          = linspace(theta_ref_r - l0_ref_degree_top, theta_ref_l + l0_ref_degree_top, round(arc_length_top/l0_ref_degree_top))';
r_top              = r_ext*ones(size(theta_top,1),1);
p_fix_ref_top_pol  = [theta_top r_top];
p_fix_ref_top      = polar2cartesian(p_fix_ref_top_pol);

% pfix left side ref region
r_length_left      = abs((r_ext - l0_ref) - (r_ext - d_ref + l0_ref));
r_left             = linspace(r_ext - l0_ref, r_ext - d_ref + l0_ref, round(r_length_left/l0_ref))';
theta_left         = theta_ref_l*ones(size(r_left,1),1);
pfix_ref_left_pol  = [theta_left r_left];
pfix_ref_left      = polar2cartesian(pfix_ref_left_pol);

%==========================================================================
% OUTPUT DATA
%==========================================================================
if strcmp(SETTINGS.mesh,'axisym') % remove interior points with x <= 0 
    pfix_vertex(pfix_vertex(:,1) <= 0,:)       = [];
    p_fix_ref_bot(p_fix_ref_bot(:,1) <= 0,:)   = [];
    pfix_ref_right(pfix_ref_right(:,1) <= 0,:) = [];
    p_fix_ref_top(p_fix_ref_top(:,1) <= 0,:)   = [];
    pfix_ref_left(pfix_ref_left(:,1) <= 0,:)   = [];
    if sum(ismember([0 r_ext],[pfix; pfix_vertex; p_fix_ref_bot; pfix_ref_right; p_fix_ref_top; pfix_ref_left],'rows'))==0
        pfix = [[0 r_ext]; pfix];
    end
end
MESH.pfix_mesh_vertex = pfix;
MESH.pfix_ref_vertex  = pfix_vertex;
MESH.pfix             = [pfix; pfix_vertex; p_fix_ref_bot; pfix_ref_right; p_fix_ref_top; pfix_ref_left];
MESH.pbnd2            = pbnd2;
GCOORD                = [MESH.pfix;MESH.pbnd1;MESH.pbnd2;MESH.pint];
EL2NOD                = delaunay(GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_POL            = cartesian2polar(GCOORD);
EL2NOD                = EL2NOD(~(sum(ismember(EL2NOD,find(abs(GCOORD_POL(:,2)-SETTINGS.r_int) < 1e-8)),2)==3),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_POL(:,2)-SETTINGS.r_int) < 1e-8);
% elements_with_3_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 3;
% remove_elements_with_3_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_3_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------
MESH.GCOORD           = GCOORD;
MESH.EL2NOD           = EL2NOD;
MESH.q                = mesh_quality(MESH.GCOORD,MESH.EL2NOD);

% SETTINGS.show_figs = 1;
% SETTINGS.save_figs = 0;
% plot_first_guess_mesh(MESH,SETTINGS)

end % END OF FUNCTION pfix_refined_region