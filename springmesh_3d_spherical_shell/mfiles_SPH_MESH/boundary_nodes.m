function [pfix,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH)
% Usage: [pfix,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Generation of boundary nodes for a 3D spherical shell mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   pfix       : [matrix]    : coord of fixed nodes
%   pbnd1      : [matrix]    : coord of boundary nodes for inner boundary
%   pbnd2      : [matrix]    : coord of boundary nodes for outer boundary
%
% JMT May 2016
% JMT Jun 2017: Fixed nodes may be created along the Z axis
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
deg2rad      = pi/180;
r_int        = SETTINGS.r_int;      % inner radius (km) of spherical shell
r_ext        = SETTINGS.r_ext;      % outer radius (km) of spherical shell
theta0       = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0         = GUIDE_MESH.phi0;     % longitude (degrees) of the point around which the refined and transition zones are defined
w_tran_deg   = GUIDE_MESH.w_tran/(deg2rad*r_ext); % width of transition zone in degrees (North-South)
theta_tran_n = theta0-w_tran_deg/2; % colatitude of the northern boundary in the transition zone
theta_tran_s = theta0+w_tran_deg/2; % colatitude of the southern boundary in the transition zone
l_tran_deg   = GUIDE_MESH.l_tran/(deg2rad*r_ext); % length of transition zone in degrees (East-West)
phi_tran_e   = phi0+l_tran_deg/2;   % longitude of the eastern boundary in the transition zone
phi_tran_w   = phi0-l_tran_deg/2;   % longitude of the western boundary in the transition zone
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees (North-South)
theta_ref_n  = theta0-w_ref_deg/2;  % colatitude of the northern boundary in the refined zone
theta_ref_s  = theta0+w_ref_deg/2;  % colatitude of the southern boundary in the refined zone
l_ref_deg    = GUIDE_MESH.l_ref/(deg2rad*r_ext);  % length of refined zone in degrees (East-West)
phi_ref_e    = phi0+l_ref_deg/2;    % longitude of the eastern boundary in the refined zone
phi_ref_w    = phi0-l_ref_deg/2;    % longitude of the western boundary in the refined zone

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

%==========================================================================
% BOUNDARY 1 (INNER BOUNDARY)
%==========================================================================
% COARSE ZONE
pbnd1_coarse         = regular_bnd_nodes(r_int,GUIDE_MESH.l0_coarse); 
[pfix1,pbnd1_coarse] = create_pfix(pbnd1_coarse,r_int);
pbnd1_coarse_sph     = cartesian2spherical(pbnd1_coarse);
theta1_coarse        = pbnd1_coarse_sph(:,1);
phi1_coarse          = pbnd1_coarse_sph(:,2);
pbnd1_coarse(theta1_coarse > theta_tran_n & theta1_coarse < theta_tran_s & ...
             phi1_coarse   > phi_tran_w   & phi1_coarse   < phi_tran_e,:) = []; % remove bnd nodes that are in transition zone
pbnd1_coarse         = reject_points_surf(pbnd1_coarse,GUIDE_MESH); % reject pbnd1_coarse using probability

% TRANSITION ZONE
pbnd1_tran           = regular_bnd_nodes(r_int,GUIDE_MESH.l0_coarse);
[~,pbnd1_tran]       = create_pfix(pbnd1_tran,r_int); % remove pfix1 values from pbnd1_tran
pbnd1_tran_sph       = cartesian2spherical(pbnd1_tran);
theta1_tran          = pbnd1_tran_sph(:,1);
phi1_tran            = pbnd1_tran_sph(:,2);
pbnd1_tran           = pbnd1_tran(theta1_tran > theta_tran_n & theta1_tran < theta_tran_s & ...
                                  phi1_tran   > phi_tran_w   & phi1_tran   < phi_tran_e,:); % take only the nodes in the transition zone
pbnd1_tran           = reject_points_surf(pbnd1_tran,GUIDE_MESH); % reject pbnd1_tran using probability

if GUIDE_MESH.d_ref == GUIDE_MESH.d_tran
    theta1_tran_temp = theta1_tran(theta1_tran > theta_tran_n & theta1_tran < theta_tran_s & ...
                                   phi1_tran   > phi_tran_w   & phi1_tran   < phi_tran_e,:); % take theta angles inside the transition zone
    phi1_tran_temp   = phi1_tran(  theta1_tran > theta_tran_n & theta1_tran < theta_tran_s & ...
                                   phi1_tran   > phi_tran_w   & phi1_tran   < phi_tran_e,:); % take phi angles inside the transition zone
    pbnd1_tran(theta1_tran_temp > theta_ref_n & theta1_tran_temp < theta_ref_s & ...
               phi1_tran_temp   > phi_ref_w   & phi1_tran_temp   < phi_ref_e,:) = []; % remove bnd nodes that are inside the refined zone
    
    % REFINED ZONE
    pbnd1_ref        = regular_bnd_nodes(r_int,GUIDE_MESH.l0_ref);
    [~,pbnd1_ref]    = create_pfix(pbnd1_ref,r_int); % remove pfix1 values from pbnd1_ref
    pbnd1_ref_sph    = cartesian2spherical(pbnd1_ref);
    theta1_ref       = pbnd1_ref_sph(:,1);
    phi1_ref         = pbnd1_ref_sph(:,2);
    pbnd1_ref        = pbnd1_ref(theta1_ref > theta_ref_n & theta1_ref < theta_ref_s & ...
                                 phi1_ref   > phi_ref_w   & phi1_ref   < phi_ref_e,:); % take only the nodes in the refined zone
    pbnd1_ref        = reject_points_surf(pbnd1_ref,GUIDE_MESH); % reject pbnd1_ref using probability
else
    pbnd1_ref = [];
end

%==========================================================================
% BOUNDARY 2 (OUTER BOUNDARY)
%==========================================================================
% COARSE ZONE
[pbnd2_coarse]       = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_coarse);
[pfix2,pbnd2_coarse] = create_pfix(pbnd2_coarse,r_ext);
[pbnd2_coarse_sph]   = cartesian2spherical(pbnd2_coarse);
theta2_coarse        = pbnd2_coarse_sph(:,1);
phi2_coarse          = pbnd2_coarse_sph(:,2);
pbnd2_coarse(theta2_coarse > theta_tran_n & theta2_coarse < theta_tran_s & ...
               phi2_coarse >   phi_tran_w &   phi2_coarse <   phi_tran_e,:) = []; % remove bnd nodes that are in the transition zone
pbnd2_coarse         = reject_points_surf(pbnd2_coarse,GUIDE_MESH); % reject pbnd2_coarse using probability

% TRANSITION ZONE
[pbnd2_tran]         = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_ref);
[~,pbnd2_tran]       = create_pfix(pbnd2_tran,r_ext); % remove pfix2 values from pbnd2_tran
[pbnd2_tran_sph]     = cartesian2spherical(pbnd2_tran);
theta2_tran          = pbnd2_tran_sph(:,1);
phi2_tran            = pbnd2_tran_sph(:,2);
pbnd2_tran           = pbnd2_tran( theta2_tran > theta_tran_n & theta2_tran < theta_tran_s & ...
                                   phi2_tran   > phi_tran_w   & phi2_tran   < phi_tran_e,:); % take only the nodes in the transition zone
theta2_tran_temp     = theta2_tran(theta2_tran > theta_tran_n & theta2_tran < theta_tran_s & ...
                                   phi2_tran   > phi_tran_w   & phi2_tran   < phi_tran_e,:); % take theta angles inside the transition zone
phi2_tran_temp       = phi2_tran(  theta2_tran > theta_tran_n & theta2_tran < theta_tran_s & ...
                                   phi2_tran   > phi_tran_w   & phi2_tran   < phi_tran_e,:); % take phi angles inside the transition zone
pbnd2_tran(theta2_tran_temp > theta_ref_n & theta2_tran_temp < theta_ref_s & ...
             phi2_tran_temp >   phi_ref_w &   phi2_tran_temp <   phi_ref_e,:) = []; % remove bnd nodes that are inside the refined zone
pbnd2_tran           = reject_points_surf(pbnd2_tran,GUIDE_MESH); % reject pbnd2_tran using probability
       
% REFINED ZONE
[pbnd2_ref]          = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_ref);
[~,pbnd2_ref]        = create_pfix(pbnd2_ref,r_ext); % remove pfix2 values from pbnd2_ref
[pbnd2_ref_sph]      = cartesian2spherical(pbnd2_ref);
theta2_ref           = pbnd2_ref_sph(:,1);
phi2_ref             = pbnd2_ref_sph(:,2);
pbnd2_ref            = pbnd2_ref(theta2_ref > theta_ref_n & theta2_ref < theta_ref_s & ...
                                 phi2_ref >   phi_ref_w &   phi2_ref <   phi_ref_e,:); % take only the nodes in the refined zone
pbnd2_ref            = reject_points_surf(pbnd2_ref,GUIDE_MESH); % reject pbnd2_ref using probability

%==========================================================================
% OUPUT DATA
%==========================================================================
pfix  = [pfix1; pfix2];
pbnd1 = [pbnd1_coarse; pbnd1_tran; pbnd1_ref];
pbnd2 = [pbnd2_coarse; pbnd2_tran; pbnd2_ref];

% % %==========================================================================
% % % CREATE pfix ON THE Z AXIS
% % %==========================================================================
% % gamma      = 1; % factor to reduce the space between points on the z axis. This is to make elements having one bar on the z axis 
% % N          = round((r_ext-r_int)/(gamma*GUIDE_MESH.l0_coarse));
% % z_pos      = linspace(r_ext,r_int,N+2);
% % pfix_z_pos = [zeros(1,size(z_pos,2)); ...
% %               zeros(1,size(z_pos,2)); ...
% %               z_pos]';
% % z_neg      = linspace(-r_int,-r_ext,N+2);
% % pfix_z_neg = [zeros(1,size(z_pos,2)); ...
% %               zeros(1,size(z_pos,2)); ...
% %               z_neg]';
% % % remove points close to the poles (they are substituted by the ones computed in pfix_z_pos and pfix_z_neg
% % pfix1(1:2,:) = [];
% % pfix2(1:2,:) = [];
% % 
% % %==========================================================================
% % % OUPUT DATA
% % %==========================================================================
% % pfix  = [pfix1; pfix2; pfix_z_pos; pfix_z_neg];
% % pbnd1 = [pbnd1_coarse; pbnd1_tran; pbnd1_ref];
% % pbnd2 = [pbnd2_coarse; pbnd2_tran; pbnd2_ref];

end % END OF FUNCTION boundary_nodes

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function GCOORD = reject_points_surf(GCOORD,GUIDE_MESH)
% Usage: GCOORD = reject_points_surf(GCOORD,GUIDE_MESH)
%
% Purpose:
%   Reject points on a surface using point density (based on Persson and
%   Strang, 2004)
%
% Input:
%   GCOORD     : [matrix]    : coordinates of points
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   GCOORD     : [matrix]    : coordinates of points
%
% JMT Jan 2017

L0_GCOORD = bar_L0_guide(GCOORD,GUIDE_MESH); % compute the desired length for each point
r0_GCOORD = sqrt(2)./(L0_GCOORD.^2);         % probability to keep the points (it is proportional to L^-2 because points are on the surface)
GCOORD    = GCOORD(rand(size(GCOORD,1),1)<r0_GCOORD./max(r0_GCOORD),:); % reject points with a probability proportional to L^-2

end % END OF FUNCTION reject_points_surf