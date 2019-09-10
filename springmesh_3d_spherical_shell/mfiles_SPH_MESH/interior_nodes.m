function pint = interior_nodes(SETTINGS,GUIDE_MESH)
% Usage: pint = interior_nodes(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Generate interior nodes (first guess for where interior nodes are)
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   pint       : [matrix]    : coordinates of interior nodes
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
deg2rad      = pi/180;
r_int        = SETTINGS.r_int;       % inner radius (km) of spherical shell
r_ext        = SETTINGS.r_ext;       % outer radius (km) of spherical shell
theta0       = GUIDE_MESH.theta0;    % colatitude (degrees) of the point around which refined and transition zones are defined
phi0         = GUIDE_MESH.phi0;      % longitude (degrees) of the point around which refined and transition zones are defined
w_tran_deg   = GUIDE_MESH.w_tran/(deg2rad*r_ext); % width of transition zone in degrees (North-South)
d_tran       = GUIDE_MESH.d_tran;    % transition zone depth (km)
theta_tran_n = theta0-w_tran_deg/2;  % colatitude of the northern boundary in the transition zone
theta_tran_s = theta0+w_tran_deg/2;  % colatitude of the southern boundary in the transition zone
l_tran_deg   = GUIDE_MESH.l_tran/(deg2rad*r_ext); % length of transition zone in degrees (East-West)
phi_tran_e   = phi0+l_tran_deg/2;    % longitude of the eastern boundary in the transition zone
phi_tran_w   = phi0-l_tran_deg/2;    % longitude of the western boundary in the transition zone
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees (North-South)
d_ref        = GUIDE_MESH.d_ref;     % refined zone depth (km)
theta_ref_n  = theta0-w_ref_deg/2;   % colatitude of the northern boundary in the refined zone
theta_ref_s  = theta0+w_ref_deg/2;   % colatitude of the southern boundary in the refined zone
l_ref_deg    = GUIDE_MESH.l_ref/(deg2rad*r_ext);  % length of refined zone in degrees (East-West)
phi_ref_e    = phi0+l_ref_deg/2;     % longitude of the eastern boundary in the refined zone
phi_ref_w    = phi0-l_ref_deg/2;     % longitude of the western boundary in the refined zone
l0_coarse    = GUIDE_MESH.l0_coarse; % desired length (km) for coarse zone
l0_tran      = GUIDE_MESH.l0_ref;    % initial guess for desired length (km) for transition zone (same as in refined zone)
l0_ref       = GUIDE_MESH.l0_ref;    % desired length (km) for refined zone

%==========================================================================
% COARSE ZONE
%==========================================================================
pint_coarse     = hexagonal_close_packing_spheres(-r_ext,r_ext,-r_ext,r_ext,-r_ext,r_ext,l0_coarse/2);
% Remove those nodes that are outside the spherical shell
r               = sqrt(pint_coarse(:,1).^2 + pint_coarse(:,2).^2 + pint_coarse(:,3).^2); % radial distance from the origin (0,0,0) to each node
height_coarse   = (sqrt(6)/3)*l0_coarse; % height of a regular tetrahedron with side l0_coarse
pint_coarse     = pint_coarse(r > r_int + height_coarse/1.5 & r < r_ext - height_coarse/1.5,:); % pint_coarse inside the spherical shell
% Remove those nodes that are inside the transition zone
pint_coarse_sph = cartesian2spherical(pint_coarse);
theta_coarse    = pint_coarse_sph(:,1);
phi_coarse      = pint_coarse_sph(:,2);
r_coarse        = pint_coarse_sph(:,3);
pint_coarse(theta_coarse > theta_tran_n & theta_coarse < theta_tran_s & ...
            phi_coarse   >   phi_tran_w & phi_coarse   <   phi_tran_e & ...
            r_coarse     > r_ext-d_tran,:) = []; % remove int nodes inside the transition zone
pint_coarse     = reject_points_vol(pint_coarse,GUIDE_MESH); % reject pint_coarse using probability

%==========================================================================
% TRANSITION ZONE
%==========================================================================
pint_tran       = hexagonal_close_packing_spheres(-r_ext,r_ext,-r_ext,r_ext,-r_ext,r_ext,l0_tran/2);
% Remove those nodes that are outside the spherical shell
r               = sqrt(pint_tran(:,1).^2 + pint_tran(:,2).^2 + pint_tran(:,3).^2); % radial distance from the origin (0,0,0) to each node
height_tran     = (sqrt(6)/3)*l0_tran; % height of a regular tetrahedron with side l0_tran
pint_tran       = pint_tran(r > r_int + height_tran/1.5 & r < r_ext - height_tran/1.5,:); % pint_tran inside the spherical shell
% Take those nodes that are inside the transition zone
pint_tran_sph   = cartesian2spherical(pint_tran);
theta_tran      = pint_tran_sph(:,1);
phi_tran        = pint_tran_sph(:,2);
r_tran          = pint_tran_sph(:,3);
pint_tran       = pint_tran(theta_tran > theta_tran_n   + height_tran/(deg2rad*(r_ext+r_int)/2) & theta_tran < theta_tran_s - height_tran/(deg2rad*(r_ext+r_int)/2) & ...
                            phi_tran   > phi_tran_w     + height_tran/(deg2rad*(r_ext+r_int)/2) & phi_tran   < phi_tran_e   - height_tran/(deg2rad*(r_ext+r_int)/2) & ...
                            r_tran     > r_ext - d_tran + height_tran/1.5                       & r_tran     < r_ext        - height_tran/1.5,:);
% Remove those nodes that are inside the refined zone
pint_tran_sph   = cartesian2spherical(pint_tran);
theta_tran_temp = pint_tran_sph(:,1);
phi_tran_temp   = pint_tran_sph(:,2);
r_tran_temp     = pint_tran_sph(:,3);
pint_tran(theta_tran_temp > theta_ref_n   - height_tran/(deg2rad*(r_ext+r_int)/2) & theta_tran_temp < theta_ref_s + height_tran/(deg2rad*(r_ext+r_int)/2) & ...
          phi_tran_temp   > phi_ref_w     - height_tran/(deg2rad*(r_ext+r_int)/2) & phi_tran_temp   < phi_ref_e   + height_tran/(deg2rad*(r_ext+r_int)/2) & ...
          r_tran_temp     > r_ext - d_ref - height_tran/1.5 ,:) = []; 
pint_tran       = reject_points_vol(pint_tran,GUIDE_MESH); % reject pint_tran using probability

%==========================================================================
% REFINED ZONE
%==========================================================================
pint_ref        = hexagonal_close_packing_spheres(-r_ext,r_ext,-r_ext,r_ext,-r_ext,r_ext,l0_ref/2);
% Remove those nodes that are outside the spherical shell
r               = sqrt(pint_ref(:,1).^2 + pint_ref(:,2).^2 + pint_ref(:,3).^2); % radial distance from the origin (0,0,0) to each node
height_ref      = (sqrt(6)/3)*l0_ref; % height of a regular tetrahedron with side l0_ref
pint_ref        = pint_ref(r > r_int + height_ref/1.5 & r < r_ext - height_ref/1.5,:); % pint_ref inside the spherical shell
% Take those nodes that are inside the refined zone
pint_ref_sph    = cartesian2spherical(pint_ref);
theta_ref       = pint_ref_sph(:,1);
phi_ref         = pint_ref_sph(:,2);
r_ref           = pint_ref_sph(:,3);
pint_ref        = pint_ref(theta_ref > theta_ref_n                    & theta_ref < theta_ref_s            & ...
                           phi_ref   > phi_ref_w                      & phi_ref   < phi_ref_e              & ...
                           r_ref     > r_ext - d_ref + height_ref/1.5 & r_ref     < r_ext - height_ref/1.5,:);
pint_ref        = reject_points_vol(pint_ref,GUIDE_MESH); % reject pint_ref using probability

%==========================================================================
% OUPUT DATA
%==========================================================================
pint = [pint_coarse; pint_tran; pint_ref];

end  % END OF SUBFUNCTION interior_nodes

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function GCOORD = reject_points_vol(GCOORD,GUIDE_MESH)
% Usage: GCOORD = reject_points_vol(GCOORD,GUIDE_MESH)
%
% Purpose:
%   Reject points inside a volume using point density (based on Persson and
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
r0_GCOORD = sqrt(2)./(L0_GCOORD.^3);         % probability to keep the points (it is proportional to L^-3 because points are in a 3D volume)
GCOORD    = GCOORD(rand(size(GCOORD,1),1)<r0_GCOORD./max(r0_GCOORD),:); % reject points with a probability proportional to L^-3

end % END OF FUNCTION reject_points_vol