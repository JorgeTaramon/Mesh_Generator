function [GCOORD] = polar2cartesian_rad(GCOORD_POL)
% Usage: [GCOORD] = polar2cartesian_rad(GCOORD_POL)
%
% Purpose: Transform polar coordinates into Cartesian coordinates.
%
% Input:
%   GCOORD_POL : [matrix] : Polar coordinates (theta,r) (in rad and km)
%
% Output:
%   GCOORD     : [matrix] : Cartesian coordinates (x,z)(in km)
%
% Colatitude (theta) is measured from +z axis in clockwise direction (range
% 0 to 2*pi). 
%
%               z
%               |     +
%               |    /
%               |th /
%               |  / r
%               |^/
%               |/
% -----------------------------> x
%               |
%               |
%               |
%               |
%               |
%               |
%
%
% JMT Apr 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

x      = GCOORD_POL(2,:).*sin(GCOORD_POL(1,:));
z      = GCOORD_POL(2,:).*cos(GCOORD_POL(1,:));

GCOORD = [x; z];

end % END OF FUNCTION polar2cartesian_rad