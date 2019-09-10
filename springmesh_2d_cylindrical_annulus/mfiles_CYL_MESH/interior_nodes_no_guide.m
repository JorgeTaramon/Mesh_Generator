function [pint] = interior_nodes_no_guide(SETTINGS)
% Usage: [pint] = interior_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generate interior nodes (first guess for where interior nodes are) for
%   a 2D cylindrical annulus regular mesh
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   pint     : [matrix]    : coordinates of interior nodes
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
r_int = SETTINGS.r_int; % inner radius (km) of cylindrical annulus
r_ext = SETTINGS.r_ext; % outer radius (km) of cylindrical annulus
h0    = SETTINGS.h0;    % desired spring length (km) for a regular mesh (no guide mesh)

%==========================================================================
% INTERIOR NODES
%==========================================================================
[pint]      = hexagonal_close_packing_circles(-r_ext,r_ext,-r_ext,r_ext,h0/2);
% Remove those nodes that are outside the cylindrical annulus
r           = sqrt(pint(:,1).^2 + pint(:,2).^2); % radial distance from the origin (0,0) to each node
pint        = pint(r > r_int + h0/4 & r < r_ext - h0/4,:); % pint inside the cylindrical annulus

if strcmp(SETTINGS.mesh,'axisym')
    pint(pint(:,1) <= 0,:) = []; % remove interior points with x <= 0 
end

end % END OF FUNCTION interior_nodes_no_guide