function [pint] = interior_nodes(SETTINGS,GUIDE_MESH)
% Usage: [pint] = interior_nodes(SETTINGS,GUIDE_MESH)
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
r_int        = SETTINGS.r_int;       % inner radius (km) of cylindrical annulus
r_ext        = SETTINGS.r_ext;       % outer radius (km) of cylindrical annulus
theta0       = GUIDE_MESH.theta0;    % colatitude (degrees) of the point around which the 
                                     % refined and transition zones are defined
w_tran_deg   = GUIDE_MESH.w_tran/(deg2rad*r_ext); % width of transition zone in degrees
d_tran       = GUIDE_MESH.d_tran;    % transition zone depth (km)
theta_tran_l = theta0-w_tran_deg/2;  % colatitude of the left boundary in the transition zone
theta_tran_r = theta0+w_tran_deg/2;  % colatitude of the right boundary in the transition zone
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees
d_ref        = GUIDE_MESH.d_ref;     % refined zone depth (km)
theta_ref_l  = theta0-w_ref_deg/2;   % colatitude of the left boundary in the refined zone
theta_ref_r  = theta0+w_ref_deg/2;   % colatitude of the right boundary in the refined zone
l0_coarse    = GUIDE_MESH.l0_coarse; % desired length (km) for coarse zone
l0_tran      = GUIDE_MESH.l0_ref;    % desired length (km) for transition zone (same as in refined zone)
l0_ref       = GUIDE_MESH.l0_ref;    % desired length (km) for refined zone


%==========================================================================
% COARSE ZONE
%==========================================================================
[pint_coarse]   = hexagonal_close_packing_circles(-r_ext,r_ext,-r_ext,r_ext,l0_coarse/2);
% Remove those nodes that are outside the cylindrical annulus
r               = sqrt(pint_coarse(:,1).^2 + pint_coarse(:,2).^2); % radial distance from the origin (0,0) to each node
pint_coarse     = pint_coarse(r > r_int + l0_coarse/4 &...
                              r < r_ext - l0_coarse/4,:); % pint_coarse inside the cylindrical annulus
% Remove those nodes that are inside the transition zone
pint_coarse_pol = cartesian2polar(pint_coarse);
theta_coarse    = pint_coarse_pol(:,1);
r_coarse        = pint_coarse_pol(:,2);
pint_coarse(r_coarse     > r_ext-d_tran & ...
            theta_coarse > theta_tran_l & theta_coarse < theta_tran_r,:) = []; % remove int nodes inside the transition zone
pint_coarse     = reject_points_surf(pint_coarse,GUIDE_MESH); % reject pint_coarse using probability

%==========================================================================
% TRANSITION ZONE
%==========================================================================
[pint_tran]     = hexagonal_close_packing_circles(-r_ext,r_ext,-r_ext,r_ext,l0_tran/2);
% Remove those nodes that are outside the cylindrical annulus
r               = sqrt(pint_tran(:,1).^2 + pint_tran(:,2).^2); % radial distance from the origin (0,0) to each node
pint_tran       = pint_tran(r > r_int + l0_tran/2 & r < r_ext - l0_tran/2,:); % pint_tran inside the cylindrical annulus
% Take those nodes that are inside the transition zone
pint_tran_pol   = cartesian2polar(pint_tran);
theta_tran      = pint_tran_pol(:,1);
r_tran          = pint_tran_pol(:,2);
pint_tran       = pint_tran(r_tran     > r_ext-d_tran + l0_tran/4                 & ...
                            r_tran     < r_ext        - l0_tran/2                 & ...
                            theta_tran > theta_tran_l + l0_tran/(4*deg2rad*r_ext) & ...
                            theta_tran < theta_tran_r - l0_tran/(4*deg2rad*r_ext),:);
% Remove those nodes that are inside the refined zone
pint_tran_pol   = cartesian2polar(pint_tran);
theta_tran_temp = pint_tran_pol(:,1);
r_tran_temp     = pint_tran_pol(:,2);
pint_tran(r_tran_temp     > r_ext-d_ref - 2*l0_ref                 & ...
          theta_tran_temp > theta_ref_l - 2*l0_ref/(r_ext*deg2rad) &...
          theta_tran_temp < theta_ref_r + 2*l0_ref/(r_ext*deg2rad),:) = []; 
pint_tran       = reject_points_surf(pint_tran,GUIDE_MESH); % reject pint_tran using probability

%==========================================================================
% REFINED ZONE
%==========================================================================
[pint_ref]      = hexagonal_close_packing_circles(-r_ext,r_ext,-r_ext,r_ext,l0_ref/2);
% Remove those nodes that are outside the cylindrical annulus
r               = sqrt(pint_ref(:,1).^2 + pint_ref(:,2).^2); % radial distance from the origin (0,0) to each node
pint_ref        = pint_ref(r > r_int + l0_ref/2 & r < r_ext - l0_ref/2,:); % pint_ref inside the cylindrical annulus
% Take those nodes that are inside the refined zone
pint_ref_pol    = cartesian2polar(pint_ref);
theta_ref       = pint_ref_pol(:,1);
r_ref           = pint_ref_pol(:,2);
pint_ref        = pint_ref(r_ref     > r_ext-d_ref + l0_ref                 & ...
                           r_ref     < r_ext       - l0_ref                 & ...
                           theta_ref > theta_ref_l + l0_ref/(r_ext*deg2rad) & ...
                           theta_ref < theta_ref_r - l0_ref/(r_ext*deg2rad),:);
pint_ref        = reject_points_surf(pint_ref,GUIDE_MESH); % reject pint_ref using probability

%==========================================================================
% OUPUT DATA
%==========================================================================
pint = [pint_coarse; pint_tran; pint_ref];

if strcmp(SETTINGS.mesh,'axisym')
    pint(pint(:,1) <= 0,:) = []; % remove interior points with x <= 0 
end

end  % END OF SUBFUNCTION interior_nodes

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