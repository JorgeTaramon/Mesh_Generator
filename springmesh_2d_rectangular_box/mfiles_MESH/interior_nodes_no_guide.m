function [pint] = interior_nodes_no_guide(SETTINGS)
% Usage: [pint] = interior_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generate interior nodes (first guess for where interior nodes are) for
%   a 2D rectangular regular mesh
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   pint     : [matrix]    : coordinates of interior nodes
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth  = SETTINGS.depth;  % domain depth (km)
length = SETTINGS.length; % domain length (km)
x0     = SETTINGS.x0;     % point around which the domain is created
z0     = SETTINGS.z0;     % point around which the domain is created
h0     = SETTINGS.h0;     % desired spring length (km) for a regular mesh

%==========================================================================
% INTERIOR NODES
%==========================================================================
xmin   = x0 - length/2;
xmax   = x0 + length/2;
zmin   = z0 - depth;
zmax   = z0;
pint   = hexagonal_close_packing_circles(xmin,xmax,zmin,zmax,h0/2);
% Remove those nodes that are very close to the boundaries
pint   = pint(pint(:,1) > xmin + h0/4 & ...
              pint(:,1) < xmax - h0/4 & ...
              pint(:,2) > zmin + h0/4 & ...
              pint(:,2) < zmax - h0/4,:);
          
end % END OF FUNCTION interior_nodes_no_guide