function [pfix,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS)
% Usage: [pfix,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generation of boundary nodes for a 3D spherical regular shell mesh.
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   pfix     : [matrix]    : coord of fixed nodes
%   pbnd1    : [matrix]    : coord of boundary nodes for inner boundary
%   pbnd2    : [matrix]    : coord of boundary nodes for outer boundary
%
% JMT May 2016
% JMT Jun 2017: Fixed nodes may be created along the Z axis
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
r_int = SETTINGS.r_int;
r_ext = SETTINGS.r_ext;
h0    = SETTINGS.h0;

%==========================================================================
% BOUNDARY 1 (INNER BOUNDARY)
%==========================================================================
if r_int > 0
    [pbnd1]       = regular_bnd_nodes(r_int,h0);
    [pfix1,pbnd1] = create_pfix(pbnd1,r_int);
elseif r_int == 0
    pfix1 = [];
    pbnd1 = [];
else
    error('r_int must be >= 0')
end

%==========================================================================
% BOUNDARY 2 (OUTER BOUNDARY)
%==========================================================================
[pbnd2]       = regular_bnd_nodes(r_ext,h0);
[pfix2,pbnd2] = create_pfix(pbnd2,r_ext);

%==========================================================================
% OUPUT DATA
%==========================================================================
pfix          = [pfix1; pfix2];

% % %==========================================================================
% % % CREATE pfix ON THE Z AXIS
% % %==========================================================================
% % gamma      = 1; % factor to reduce the space between points on the z axis. This is to make elements having one bar on the z axis 
% % N          = round((r_ext-r_int)/(gamma*h0));
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
% % pfix          = [pfix1; pfix2; pfix_z_pos; pfix_z_neg];

end % END OF FUNCTION boundary_nodes_no_guide