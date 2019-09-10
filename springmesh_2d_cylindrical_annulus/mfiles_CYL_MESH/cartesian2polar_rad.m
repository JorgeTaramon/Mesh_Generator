function [GCOORD_POL] = cartesian2polar_rad(GCOORD)
% Usage: [GCOORD_POL] = cartesian2polar_rad(GCOORD)
%
% Purpose: Transform Cartesian coordinates into Polar coordinates.
%
% Input:
%   GCOORD     : [matrix] : Cartesian coordinates (x,z)(in km)
%
% Output:
%   GCOORD_POL : [matrix] : Polar coordinates (theta,r) (in rad and km)
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
% JMT Feb 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

r              = sqrt(GCOORD(1,:).^2+GCOORD(2,:).^2);
theta          = atan2(GCOORD(1,:),GCOORD(2,:));
theta(theta<0) = theta(theta<0)+2*pi;

GCOORD_POL     = [theta; r];

end % END OF FUNCTION cartesian2polar_rad