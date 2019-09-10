function MESH = fixed_nodes_refined_region(MESH,GUIDE_MESH,SETTINGS)
% Usage: MESH = fixed_nodes_refined_region(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Define pbnd for the embedded high resolution sub-region. These nodes
%   will move parallel to the boundaries where they have been defined.
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
r_ext        = SETTINGS.r_ext;      % outer radius (km) of spherical shell
theta0       = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0         = GUIDE_MESH.phi0;     % longitude (degrees) of the point around which the refined and transition zones are defined
d_ref        = GUIDE_MESH.d_ref;    % refined zone depth (km)
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees (North-South)
theta_ref_n  = theta0-w_ref_deg/2;  % colatitude of the northern boundary in the refined zone
theta_ref_s  = theta0+w_ref_deg/2;  % colatitude of the southern boundary in the refined zone
l_ref_deg    = GUIDE_MESH.l_ref/(deg2rad*r_ext);  % length of refined zone in degrees (East-West)
phi_ref_e    = phi0+l_ref_deg/2;    % longitude of the eastern boundary in the refined zone
phi_ref_w    = phi0-l_ref_deg/2;    % longitude of the western boundary in the refined zone
l0_ref       = GUIDE_MESH.l0_ref;   % desired spring length (km) for refined zone
pfix         = MESH.pfix;
pbnd2        = MESH.pbnd2;

% Top  view of transition and refined zones
%
%
%   l0_coarse                                       l0_coarse
%          ._________________________________________. - - - - -> theta_tran_n
%          |                                         |
%          |                                         |
%          |             TRANSITION ZONE             |
%          |                                         |
%          |                                         |
%          |      l0_ref._______________.l0_ref      | - - - - -> theta_ref_n
%          |            |               |            |
%          |            |  REFINED ZONE |            |
%          |            |               |            |
%          |            |       .       |            |
%          |            | (theta0,phi0) |            |
%          |            |               |            |
%          |      l0_ref._______________.l0_ref      | - - - - -> theta_ref_s
%          |                                         |
%          |                                         |
%          |                                         |
%          |                                         |
%          |                                         |
%          ._________________________________________. - - - - -> theta_tran_s
%   l0_coarse                                       l0_coarse
%          |            |               |            |
%
%          |            |               |            |
%
%          |            |               |            |
%   phi_tran_w       phi_ref_w       phi_ref_e    phi_tran_e
%

%===============================================================================
% REMOVE PFIX AND PBND ON BOUNDARY 2 (OUTER BOUNDARY) INSIDE THE REFINED REGION
%===============================================================================
pfix_sph          = cartesian2spherical(pfix);
pfix(pfix_sph(:,1) >= theta_ref_n & pfix_sph(:,1) <= theta_ref_s & ...
     pfix_sph(:,2) >= phi_ref_w   & pfix_sph(:,2) <= phi_ref_e   & ...
     pfix_sph(:,3) > r_ext - 1e-8,:) = [];

l0_ref_degree_top = l0_ref/(r_ext*deg2rad);
pbnd2_sph         = cartesian2spherical(pbnd2);
pbnd2(pbnd2_sph(:,1) > theta_ref_n - l0_ref_degree_top & pbnd2_sph(:,1) < theta_ref_s + l0_ref_degree_top & ...
      pbnd2_sph(:,2) > phi_ref_w - l0_ref_degree_top   & pbnd2_sph(:,2) < phi_ref_e + l0_ref_degree_top   & ...
      pbnd2_sph(:,3) > r_ext - 1e-8,:) = [];

%==========================================================================
% DEFINE PFIX VERTICES REFINED REGION
%==========================================================================
pfix_vertex_sph    = [theta_ref_s phi_ref_w  r_ext - d_ref ; ...
                      theta_ref_s phi_ref_e  r_ext - d_ref ; ...
                      theta_ref_n phi_ref_e  r_ext - d_ref ; ...
                      theta_ref_n phi_ref_w  r_ext - d_ref ; ...
                      theta_ref_s phi_ref_w  r_ext ; ...
                      theta_ref_s phi_ref_e  r_ext ; ...
                      theta_ref_n phi_ref_e  r_ext ; ...
                      theta_ref_n phi_ref_w  r_ext ];
pfix_vertex        = spherical2cartesian(pfix_vertex_sph);

%==========================================================================
% DEFINE PFIX EDGES REFINED REGION (12 EDGES)
%==========================================================================

% BOTTOM EDGES (4 EDGES)
l0_ref_degree_bot    = l0_ref/((r_ext - d_ref)*deg2rad);
% edge 1 (from (theta_ref_s,phi_ref_w,r_ref_bot) to (theta_ref_s,phi_ref_e,r_ref_bot))
arc_length_bot_1     = abs((phi_ref_w + l0_ref_degree_bot) - (phi_ref_e - l0_ref_degree_bot));
phi_bot_1            = linspace(phi_ref_w + l0_ref_degree_bot, phi_ref_e - l0_ref_degree_bot, round(arc_length_bot_1/l0_ref_degree_bot))';
r_bot_1              = (r_ext - d_ref)*ones(size(phi_bot_1,1),1);
theta_bot_1          = theta_ref_s*ones(size(phi_bot_1,1),1);
pfix_ref_bot_1_sph   = [theta_bot_1 phi_bot_1 r_bot_1];
pfix_ref_bot_1       = spherical2cartesian(pfix_ref_bot_1_sph);
% edge 2 (from (theta_ref_s,phi_ref_e,r_ref_bot) to (theta_ref_n,phi_ref_e,r_ref_bot))
arc_length_bot_2     = abs((theta_ref_s - l0_ref_degree_bot) - (theta_ref_n + l0_ref_degree_bot));
theta_bot_2          = linspace(theta_ref_s - l0_ref_degree_bot, theta_ref_n + l0_ref_degree_bot, round(arc_length_bot_2/l0_ref_degree_bot))';
phi_bot_2            = phi_ref_e*ones(size(theta_bot_2,1),1);
r_bot_2              = (r_ext - d_ref)*ones(size(theta_bot_2,1),1);
pfix_ref_bot_2_sph   = [theta_bot_2 phi_bot_2 r_bot_2];
pfix_ref_bot_2       = spherical2cartesian(pfix_ref_bot_2_sph);
% edge 3 (from (theta_ref_n,phi_ref_e,r_ref_bot) to (theta_ref_n,phi_ref_w,r_ref_bot))
arc_length_bot_3     = abs((phi_ref_e - l0_ref_degree_bot) - (phi_ref_w + l0_ref_degree_bot));
phi_bot_3            = linspace(phi_ref_e - l0_ref_degree_bot, phi_ref_w + l0_ref_degree_bot, round(arc_length_bot_3/l0_ref_degree_bot))';
r_bot_3              = (r_ext - d_ref)*ones(size(phi_bot_3,1),1);
theta_bot_3          = theta_ref_n*ones(size(phi_bot_3,1),1);
pfix_ref_bot_3_sph   = [theta_bot_3 phi_bot_3 r_bot_3];
pfix_ref_bot_3       = spherical2cartesian(pfix_ref_bot_3_sph);
% edge 4 (from (theta_ref_n,phi_ref_w,r_ref_bot) to (theta_ref_s,phi_ref_w,r_ref_bot))
arc_length_bot_4     = abs((theta_ref_n + l0_ref_degree_bot) - (theta_ref_s - l0_ref_degree_bot));
theta_bot_4          = linspace(theta_ref_n + l0_ref_degree_bot, theta_ref_s - l0_ref_degree_bot, round(arc_length_bot_4/l0_ref_degree_bot))';
phi_bot_4            = phi_ref_w*ones(size(theta_bot_4,1),1);
r_bot_4              = (r_ext - d_ref)*ones(size(theta_bot_4,1),1);
pfix_ref_bot_4_sph   = [theta_bot_4 phi_bot_4 r_bot_4];
pfix_ref_bot_4       = spherical2cartesian(pfix_ref_bot_4_sph);

% TOP EDGES (4 EDGES)
% edge 1 (from (theta_ref_s,phi_ref_w,r_ext) to (theta_ref_s,phi_ref_e,r_ext))
arc_length_top_1     = abs((phi_ref_w + l0_ref_degree_top) - (phi_ref_e - l0_ref_degree_top));
phi_top_1            = linspace(phi_ref_w + l0_ref_degree_top, phi_ref_e - l0_ref_degree_top, round(arc_length_top_1/l0_ref_degree_top))';
r_top_1              = r_ext*ones(size(phi_top_1,1),1);
theta_top_1          = theta_ref_s*ones(size(phi_top_1,1),1);
pfix_ref_top_1_sph   = [theta_top_1 phi_top_1 r_top_1];
pfix_ref_top_1       = spherical2cartesian(pfix_ref_top_1_sph);
% edge 2 (from (theta_ref_s,phi_ref_e,r_ext) to (theta_ref_n,phi_ref_e,r_ext))
arc_length_top_2     = abs((theta_ref_s - l0_ref_degree_top) - (theta_ref_n + l0_ref_degree_top));
theta_top_2          = linspace(theta_ref_s - l0_ref_degree_top, theta_ref_n + l0_ref_degree_top, round(arc_length_top_2/l0_ref_degree_top))';
phi_top_2            = phi_ref_e*ones(size(theta_top_2,1),1);
r_top_2              = r_ext*ones(size(theta_top_2,1),1);
pfix_ref_top_2_sph   = [theta_top_2 phi_top_2 r_top_2];
pfix_ref_top_2       = spherical2cartesian(pfix_ref_top_2_sph);
% edge 3 (from (theta_ref_n,phi_ref_e,r_ext) to (theta_ref_n,phi_ref_w,r_ext))
arc_length_top_3     = abs((phi_ref_e - l0_ref_degree_top) - (phi_ref_w + l0_ref_degree_top));
phi_top_3            = linspace(phi_ref_e - l0_ref_degree_top, phi_ref_w + l0_ref_degree_top, round(arc_length_top_3/l0_ref_degree_top))';
r_top_3              = r_ext*ones(size(phi_top_3,1),1);
theta_top_3          = theta_ref_n*ones(size(phi_top_3,1),1);
pfix_ref_top_3_sph   = [theta_top_3 phi_top_3 r_top_3];
pfix_ref_top_3       = spherical2cartesian(pfix_ref_top_3_sph);
% edge 4 (from (theta_ref_n,phi_ref_w,r_ext) to (theta_ref_s,phi_ref_w,r_ext))
arc_length_top_4     = abs((theta_ref_n + l0_ref_degree_top) - (theta_ref_s - l0_ref_degree_top));
theta_top_4          = linspace(theta_ref_n + l0_ref_degree_top, theta_ref_s - l0_ref_degree_top, round(arc_length_top_4/l0_ref_degree_top))';
phi_top_4            = phi_ref_w*ones(size(theta_top_4,1),1);
r_top_4              = r_ext*ones(size(theta_top_4,1),1);
pfix_ref_top_4_sph   = [theta_top_4 phi_top_4 r_top_4];
pfix_ref_top_4       = spherical2cartesian(pfix_ref_top_4_sph);

% RADIAL EDGES
r_length             = abs((r_ext - d_ref + l0_ref/2) - (r_ext - l0_ref/2));
r                    = linspace(r_ext - d_ref + l0_ref, r_ext - l0_ref, round(r_length/l0_ref))';
% edge 1 (from (theta_ref_s,phi_ref_w,r_ref_bot) to (theta_ref_s,phi_ref_w,r_ext))
theta_1              = theta_ref_s*ones(size(r,1),1);
phi_1                = phi_ref_w*ones(size(r,1),1);
pfix_ref_edge_1_sph  = [theta_1 phi_1 r];
pfix_ref_edge_1      = spherical2cartesian(pfix_ref_edge_1_sph);
% edge 2 (from (theta_ref_s,phi_ref_e,r_ref_bot) to (theta_ref_s,phi_ref_e,r_ext))
theta_2              = theta_ref_s*ones(size(r,1),1);
phi_2                = phi_ref_e*ones(size(r,1),1);
pfix_ref_edge_2_sph  = [theta_2 phi_2 r];
pfix_ref_edge_2      = spherical2cartesian(pfix_ref_edge_2_sph);
% edge 3 (from (theta_ref_n,phi_ref_e,r_ref_bot) to (theta_ref_n,phi_ref_e,r_ext))
theta_3              = theta_ref_n*ones(size(r,1),1);
phi_3                = phi_ref_e*ones(size(r,1),1);
pfix_ref_edge_3_sph  = [theta_3 phi_3 r];
pfix_ref_edge_3      = spherical2cartesian(pfix_ref_edge_3_sph);
% edge 4 (from (theta_ref_n,phi_ref_w,r_ref_bot) to (theta_ref_n,phi_ref_w,r_ext))
theta_4              = theta_ref_n*ones(size(r,1),1);
phi_4                = phi_ref_w*ones(size(r,1),1);
pfix_ref_edge_4_sph  = [theta_4 phi_4 r];
pfix_ref_edge_4      = spherical2cartesian(pfix_ref_edge_4_sph);

%==========================================================================
% DEFINE PBND FACES REFINED REGION (6 FACES)
%==========================================================================

% BOTTOM FACE
pfix_ref_face_bot     = regular_bnd_nodes(r_ext - d_ref,l0_ref);
pfix_ref_face_bot_sph = cartesian2spherical(pfix_ref_face_bot);
theta                 = pfix_ref_face_bot_sph(:,1);
phi                   = pfix_ref_face_bot_sph(:,2);
pfix_ref_face_bot     = pfix_ref_face_bot(theta >= theta_ref_n + l0_ref_degree_bot & theta <= theta_ref_s - l0_ref_degree_bot & ...
                                          phi   >= phi_ref_w + l0_ref_degree_bot   & phi <= phi_ref_e - l0_ref_degree_bot     ,:);

% TOP FACE
pfix_ref_face_top     = regular_bnd_nodes(r_ext,l0_ref);
pfix_ref_face_top_sph = cartesian2spherical(pfix_ref_face_top);
theta                 = pfix_ref_face_top_sph(:,1);
phi                   = pfix_ref_face_top_sph(:,2);
pfix_ref_face_top     = pfix_ref_face_top(theta >= theta_ref_n + l0_ref_degree_top & theta <= theta_ref_s - l0_ref_degree_top & ...
                                          phi   >= phi_ref_w + l0_ref_degree_top   & phi <= phi_ref_e - l0_ref_degree_top     ,:);
% FACE 1
% (theta_ref_s,phi_ref_w,r_ref_bot) -> (theta_ref_s,phi_ref_w,r_ext) -> (theta_ref_s,phi_ref_e,r_ext) -> (theta_ref_s,phi_ref_e,r_ref_bot)
pfix_face_1_sph = hexagonal_close_packing_circles(r_ext*cosd(phi_ref_e - (l0_ref_degree_bot + l0_ref_degree_top)/2), ...
                                                  r_ext*cosd(phi_ref_w + (l0_ref_degree_bot + l0_ref_degree_top)/2), ...
                                                  r_ext - d_ref + l0_ref/2, ...
                                                  r_ext - l0_ref/2,          ...
                                                  l0_ref/2);
phi_face_1      = atan2d(pfix_face_1_sph(:,2),pfix_face_1_sph(:,1));
r_face_1        = pfix_face_1_sph(:,2);
% traslate coord to the center of the face
mean_phi_face_1 = (phi_ref_w + phi_ref_e)/2;
mean_r_face_1   = (r_ext - d_ref + r_ext)/2;
mean_phi        = mean(phi_face_1);
mean_r          = mean(r_face_1);
vector_phi_face_1    = mean_phi_face_1 - mean_phi;
vector_r_face_1      = mean_r_face_1 - mean_r;
phi_face_1_traslated = phi_face_1 + vector_phi_face_1;
r_face_1_traslated   = r_face_1 + vector_r_face_1;
theta_face_1    = theta_ref_s*ones(size(phi_face_1,1),1);
pfix_face_1_sph = [theta_face_1 phi_face_1_traslated r_face_1_traslated];
pfix_face_1     = spherical2cartesian(pfix_face_1_sph);

% FACE 2
% (theta_ref_s,phi_ref_e,r_ref_bot) -> (theta_ref_s,phi_ref_e,r_ext) -> (theta_ref_n,phi_ref_e,r_ext) -> (theta_ref_n,phi_ref_e,r_ref_bot)
pfix_face_2_sph = hexagonal_close_packing_circles(r_ext - d_ref + l0_ref/2, ...
                                                  r_ext - l0_ref/2,          ...
                                                  r_ext*cosd(theta_ref_s - (l0_ref_degree_bot + l0_ref_degree_top)/4), ...
                                                  r_ext*cosd(theta_ref_n + (l0_ref_degree_bot + l0_ref_degree_top)/4), ...
                                                  l0_ref/2);
theta_face_2    = atan2d(pfix_face_2_sph(:,1),pfix_face_2_sph(:,2));
r_face_2        = pfix_face_2_sph(:,1);
% traslate coord to the center of the face
mean_theta_face_2 = (theta_ref_n + theta_ref_s)/2;
mean_r_face_2   = (r_ext - d_ref + r_ext)/2;
mean_theta      = mean(theta_face_2);
mean_r          = mean(r_face_2);
vector_theta_face_2    = mean_theta_face_2 - mean_theta;
vector_r_face_2        = mean_r_face_2 - mean_r;
theta_face_2_traslated = theta_face_2 + vector_theta_face_2;
r_face_2_traslated     = r_face_2 + vector_r_face_2;
phi_face_2      = phi_ref_e*ones(size(theta_face_2,1),1);
pfix_face_2_sph = [theta_face_2_traslated phi_face_2 r_face_2_traslated];
pfix_face_2     = spherical2cartesian(pfix_face_2_sph);

% FACE 3
% (theta_ref_n,phi_ref_e,r_ref_bot) -> (theta_ref_n,phi_ref_e,r_ext) -> (theta_ref_n,phi_ref_w,r_ext) -> (theta_ref_n,phi_ref_w,r_ref_bot)
pfix_face_3_sph = hexagonal_close_packing_circles(r_ext*cosd(phi_ref_e - (l0_ref_degree_bot + l0_ref_degree_top)/2), ...
                                                  r_ext*cosd(phi_ref_w + (l0_ref_degree_bot + l0_ref_degree_top)/2), ...
                                                  r_ext - d_ref + l0_ref/2, ...
                                                  r_ext - l0_ref/2,          ...
                                                  l0_ref/2);
phi_face_3      = atan2d(pfix_face_3_sph(:,2),pfix_face_3_sph(:,1));
r_face_3        = pfix_face_3_sph(:,2);
% traslate coord to the center of the face
mean_phi_face_3 = (phi_ref_w + phi_ref_e)/2;
mean_r_face_3   = (r_ext - d_ref + r_ext)/2;
mean_phi        = mean(phi_face_3);
mean_r          = mean(r_face_3);
vector_phi_face_3    = mean_phi_face_3 - mean_phi;
vector_r_face_3      = mean_r_face_3 - mean_r;
phi_face_3_traslated = phi_face_3 + vector_phi_face_3;
r_face_3_traslated   = r_face_3 + vector_r_face_3;
theta_face_3    = theta_ref_n*ones(size(phi_face_3,1),1);
pfix_face_3_sph = [theta_face_3 phi_face_3_traslated r_face_3_traslated];
pfix_face_3     = spherical2cartesian(pfix_face_3_sph);

% FACE 4
% (theta_ref_n,phi_ref_w,r_ref_bot) -> (theta_ref_n,phi_ref_w,r_ext) -> (theta_ref_s,phi_ref_w,r_ext) -> (theta_ref_s,phi_ref_w,r_ref_bot)
pfix_face_4_sph = hexagonal_close_packing_circles(r_ext - d_ref + l0_ref/2, ...
                                                  r_ext - l0_ref/2,          ...
                                                  r_ext*cosd(theta_ref_s - (l0_ref_degree_bot + l0_ref_degree_top)/4), ...
                                                  r_ext*cosd(theta_ref_n + (l0_ref_degree_bot + l0_ref_degree_top)/4), ...
                                                  l0_ref/2);
theta_face_4    = atan2d(pfix_face_4_sph(:,1),pfix_face_4_sph(:,2));
r_face_4        = pfix_face_4_sph(:,1);
% traslate coord to the center of the face
mean_theta_face_4 = (theta_ref_n + theta_ref_s)/2;
mean_r_face_4   = (r_ext - d_ref + r_ext)/2;
mean_theta      = mean(theta_face_4);
mean_r          = mean(r_face_4);
vector_theta_face_4    = mean_theta_face_4 - mean_theta;
vector_r_face_4        = mean_r_face_4 - mean_r;
theta_face_4_traslated = theta_face_4 + vector_theta_face_4;
r_face_4_traslated     = r_face_4 + vector_r_face_4;
phi_face_4      = phi_ref_w*ones(size(theta_face_4,1),1);
pfix_face_4_sph = [theta_face_4_traslated phi_face_4 r_face_4_traslated];
pfix_face_4     = spherical2cartesian(pfix_face_4_sph);

%==========================================================================
% OUTPUT DATA
%==========================================================================
MESH.pfix_mesh_vertex  = pfix;
MESH.pfix_ref_vertex   = pfix_vertex;
MESH.pfix              = [pfix; ...
                          pfix_vertex; ...
                          pfix_ref_bot_1;  pfix_ref_bot_2;  pfix_ref_bot_3;  pfix_ref_bot_4;  ...
                          pfix_ref_top_1;  pfix_ref_top_2;  pfix_ref_top_3;  pfix_ref_top_4;  ...
                          pfix_ref_edge_1; pfix_ref_edge_2; pfix_ref_edge_3; pfix_ref_edge_4; ...
                          pfix_ref_face_bot; pfix_face_1; pfix_face_2; pfix_face_3; pfix_face_4];
MESH.pbnd2             = [pbnd2; pfix_ref_face_top]; % add pfix_ref_face_top to pbnd2 to allow them moving on the surface
GCOORD                 = [MESH.pfix;              ...
                          MESH.pbnd1;             ...
                          MESH.pbnd2;             ...
                          MESH.pint];
EL2NOD                 = delaunay(GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH             = cartesian2spherical(GCOORD);
EL2NOD                 = EL2NOD(~(sum(ismember(EL2NOD,find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8);
% elements_with_4_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 4;
% remove_elements_with_4_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_4_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------
MESH.GCOORD            = GCOORD;
MESH.EL2NOD            = EL2NOD;
MESH.q                 = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
MESH.s                 = shape_measure(MESH.GCOORD,MESH.EL2NOD);

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
    MESH.iter      = 1;
end

% figure(55)
% clf
% scatter3(pfix_vertex(:,1),pfix_vertex(:,2),pfix_vertex(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% hold on
% view(142.5,30)
% xlabel('X (km)')
% ylabel('Y (km)')
% zlabel('Z (km)');
% axis equal
% scatter3(pfix_ref_bot_1(:,1),pfix_ref_bot_1(:,2),pfix_ref_bot_1(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_bot_2(:,1),pfix_ref_bot_2(:,2),pfix_ref_bot_2(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_bot_3(:,1),pfix_ref_bot_3(:,2),pfix_ref_bot_3(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_bot_4(:,1),pfix_ref_bot_4(:,2),pfix_ref_bot_4(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% 
% scatter3(pfix_ref_top_1(:,1),pfix_ref_top_1(:,2),pfix_ref_top_1(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_top_2(:,1),pfix_ref_top_2(:,2),pfix_ref_top_2(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_top_3(:,1),pfix_ref_top_3(:,2),pfix_ref_top_3(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_top_4(:,1),pfix_ref_top_4(:,2),pfix_ref_top_4(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% 
% scatter3(pfix_ref_edge_1(:,1),pfix_ref_edge_1(:,2),pfix_ref_edge_1(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_edge_2(:,1),pfix_ref_edge_2(:,2),pfix_ref_edge_2(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_edge_3(:,1),pfix_ref_edge_3(:,2),pfix_ref_edge_3(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_edge_4(:,1),pfix_ref_edge_4(:,2),pfix_ref_edge_4(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% 
% scatter3(pfix_face_1(:,1),pfix_face_1(:,2),pfix_face_1(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% scatter3(pfix_face_2(:,1),pfix_face_2(:,2),pfix_face_2(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% scatter3(pfix_face_3(:,1),pfix_face_3(:,2),pfix_face_3(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% scatter3(pfix_face_4(:,1),pfix_face_4(:,2),pfix_face_4(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% 
% scatter3(pfix_ref_face_bot(:,1),pfix_ref_face_bot(:,2),pfix_ref_face_bot(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
% scatter3(pfix_ref_face_top(:,1),pfix_ref_face_top(:,2),pfix_ref_face_top(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])

end % END OF FUNCTION fixed_nodes_refined_region