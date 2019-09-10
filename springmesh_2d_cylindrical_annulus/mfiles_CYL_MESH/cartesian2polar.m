function [GCOORD_POL] = cartesian2polar(GCOORD)
% Usage: [GCOORD_POL] = cartesian2polar(GCOORD)
%
% Purpose: Transform Cartesian coordinates into Polar coordinates.
%
% Input:
%   GCOORD     : [matrix] : Cartesian coordinates (x,y)(in km)
%
% Output:
%   GCOORD_POL : [matrix] : Polar coordinates (theta,r) (in degrees and km)
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

r              = sqrt(GCOORD(:,1).^2+GCOORD(:,2).^2);
theta          = atan2d(GCOORD(:,1),GCOORD(:,2));
theta(theta<0) = theta(theta<0)+360;

GCOORD_POL     = [theta r];

end % END OF FUNCTION cartesian2polar