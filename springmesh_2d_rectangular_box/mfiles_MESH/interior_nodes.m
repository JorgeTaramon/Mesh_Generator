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
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth      = SETTINGS.depth;         % domain depth (km)
length     = SETTINGS.length;        % domain length (km)
x0         = GUIDE_MESH.x0;          % point around which the refined and transition zones are defined
z0         = GUIDE_MESH.z0;          % point around which the refined and transition zones are defined
d_tran     = GUIDE_MESH.d_tran;      % transition zone depth (km)
l_tran     = GUIDE_MESH.l_tran;      % length of transition zone (km)
x_tran_l   = x0 - l_tran/2;          % left boundary of the transition zone
x_tran_r   = x0 + l_tran/2;          % right boundary of the transition zone
d_ref      = GUIDE_MESH.d_ref;       % refined zone depth (km)
l_ref      = GUIDE_MESH.l_ref;       % length of refined zone (km)
x_ref_l    = x0 - l_ref/2;           % left boundary of the refined zone
x_ref_r    = x0 + l_ref/2;           % right boundary of the refined zone
l0_coarse  = GUIDE_MESH.l0_coarse;   % desired length (km) for coarse zone
l0_tran    = GUIDE_MESH.l0_ref;      % desired length (km) for transition zone (same as in refined zone)
l0_ref     = GUIDE_MESH.l0_ref;      % desired length (km) for refined zone

%==========================================================================
% COARSE ZONE
%==========================================================================
xmin_coarse   = x0 - length/2;
xmax_coarse   = x0 + length/2;
zmin_coarse   = z0 - depth;
zmax_coarse   = z0;
pint_coarse   = hexagonal_close_packing_circles(xmin_coarse,xmax_coarse,zmin_coarse,zmax_coarse,l0_coarse/2);
% Remove those nodes that are very close to the boundaries
x_coarse      = pint_coarse(:,1);
z_coarse      = pint_coarse(:,2);
pint_coarse   = pint_coarse(x_coarse > xmin_coarse + l0_coarse/4 & ...
                            x_coarse < xmax_coarse - l0_coarse/4 & ...
                            z_coarse > zmin_coarse + l0_coarse/4 & ...
                            z_coarse < zmax_coarse - l0_coarse/4,:);
% Remove those nodes that are inside the transition zone
x_coarse      = pint_coarse(:,1);
z_coarse      = pint_coarse(:,2);
pint_coarse(x_coarse > x_tran_l & ...
            x_coarse < x_tran_r & ...
            z_coarse > z0 - d_tran,:) = [];
pint_coarse   = reject_points_surf(pint_coarse,GUIDE_MESH); % reject pint_coarse using probability

%==========================================================================
% TRANSITION ZONE
%==========================================================================
xmin_tran     = x_tran_l;
xmax_tran     = x_tran_r;
zmin_tran     = z0 - d_tran;
zmax_tran     = z0;
pint_tran     = hexagonal_close_packing_circles(xmin_tran,xmax_tran,zmin_tran,zmax_tran,l0_tran/2);
% Remove those nodes that are very close to the boundaries
x_tran        = pint_tran(:,1);
z_tran        = pint_tran(:,2);
pint_tran     = pint_tran(x_tran > xmin_tran + l0_tran/2 & ...
                          x_tran < xmax_tran - l0_tran/2 & ...
                          z_tran > zmin_tran + l0_tran/2 & ...
                          z_tran < zmax_tran - l0_tran/2,:);
% Remove those nodes that are inside the refined zone
x_tran        = pint_tran(:,1);
z_tran        = pint_tran(:,2);
pint_tran(x_tran > x_ref_l - 2*l0_ref & ...
          x_tran < x_ref_r + 2*l0_ref & ...
          z_tran > z0 - d_ref - 2*l0_ref,:) = [];
pint_tran     = reject_points_surf(pint_tran,GUIDE_MESH); % reject pint_tran using probability

%==========================================================================
% REFINED ZONE
%==========================================================================
xmin_ref      = x_ref_l;
xmax_ref      = x_ref_r;
zmin_ref      = z0 - d_ref;
zmax_ref      = z0;
pint_ref      = hexagonal_close_packing_circles(xmin_ref,xmax_ref,zmin_ref,zmax_ref,l0_ref/2);
% Remove those nodes that are very close to the boundaries
x_ref         = pint_ref(:,1);
z_ref         = pint_ref(:,2);
pint_ref      = pint_ref(x_ref > xmin_ref + l0_ref & ...
                         x_ref < xmax_ref - l0_ref & ...
                         z_ref > zmin_ref + l0_ref & ...
                         z_ref < zmax_ref - l0_ref,:);
pint_ref      = reject_points_surf(pint_ref,GUIDE_MESH); % reject pint_ref using probability

%==========================================================================
% OUPUT DATA
%==========================================================================
pint = [pint_coarse; pint_tran; pint_ref];

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
% JMT Jun 2017

if isempty(GCOORD)
    return
end
L0_GCOORD = bar_L0_guide(GCOORD,GUIDE_MESH); % compute the desired length for each point
r0_GCOORD = sqrt(2)./(L0_GCOORD.^2);         % probability to keep the points (it is proportional to L^-2 because points are on the surface)
GCOORD    = GCOORD(rand(size(GCOORD,1),1)<r0_GCOORD./max(r0_GCOORD),:); % reject points with a probability proportional to L^-2

end % END OF FUNCTION reject_points_surf