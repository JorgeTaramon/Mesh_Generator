function [GCOORD] = polar2cartesian(GCOORD_POL)
% Usage: [GCOORD] = polar2cartesian(GCOORD_POL)
%
% Purpose: Transform polar coordinates into Cartesian coordinates.
%
% Input:
%   GCOORD_POL : [matrix] : Polar coordinates (theta,r) (in degrees and km)
%
% Output:
%   GCOORD     : [matrix] : Cartesian coordinates (x,y)(in km)
%
% Colatitude (theta) is measured from +y axis in clockwise direction (range
% 0 to 360). 
%
%               y
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
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

x      = GCOORD_POL(:,2).*sind(GCOORD_POL(:,1));
y      = GCOORD_POL(:,2).*cosd(GCOORD_POL(:,1));

GCOORD = [x y];

end % END OF FUNCTION polar2cartesian