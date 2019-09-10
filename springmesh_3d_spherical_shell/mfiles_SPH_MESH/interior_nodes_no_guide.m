function [pint] = interior_nodes_no_guide(SETTINGS)
% Usage: [pint] = interior_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generate interior nodes (first guess for where interior nodes are) for
%   a 3D spherical regular shell mesh
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
r_int = SETTINGS.r_int; % inner radius (km) of spherical shell
r_ext = SETTINGS.r_ext; % outer radius (km) of spherical shell
h0    = SETTINGS.h0;    % desired spring length (km) for a regular mesh (no guide mesh)

%==========================================================================
% INTERIOR NODES
%==========================================================================
[pint]      = hexagonal_close_packing_spheres(-r_ext,r_ext,-r_ext,r_ext,-r_ext,r_ext,h0/2);
r           = sqrt(pint(:,1).^2 + pint(:,2).^2 + pint(:,3).^2); % radial distance from the origin (0,0,0) to each node
if r_int > 0
    % Remove those nodes that are outside the spherical shell
    pint    = pint(r > r_int + h0/4 & r < r_ext - h0/4,:); % pint_coarse inside the spherical shell
else % means that r_int = 0
    % Remove those nodes that are outside the sphere
    pint    = pint(r < r_ext - h0/4,:); % pint_coarse inside the sphere
end

end % END OF FUNCTION interior_nodes_no_guide