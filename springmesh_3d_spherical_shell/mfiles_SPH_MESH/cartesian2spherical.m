function [GCOORD_SPH] = cartesian2spherical(GCOORD)
% Usage: [GCOORD_SPH] = cartesian2spherical(GCOORD)
%
% Purpose:
%   Transform Cartesian coordinates into spherical coordinates
%   
% Input:
%   GCOORD     : [matrix] : Cartesian coordinates (x,y,z)(in km)
%
% Output:
%   GCOORD_SPH : [matrix] : Spherical coordinates (theta,phi,r) 
%                           (in degrees,degrees and km)
%
% Colatitude (theta) is measured from +z axis to -z (range 0 to 180).
% Longitude (phi) is measured from +X axis in counterclowise direction
% (west to east) (range 0 to 360).
%
%                 Z
%                 |
%                 |_
%                 | \
%                 |  \
%                 |   \
%                 |    \
%                 |     +
%                 |    /.
%                 |th / .
%                 |  /r .
%                 |^/   .
%                 |/    .
%                _-_----.-------_----> Y
%              _-\_/\   .     _-
%            _-  phi \  .   _-
%          _-         \ . _-
%        _-____________\.-
%      _-                
%    X              
%
% JMT Jun 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

r            = sqrt(GCOORD(:,1).^2+GCOORD(:,2).^2+GCOORD(:,3).^2);
theta        = atan2d(sqrt(GCOORD(:,1).^2+GCOORD(:,2).^2),GCOORD(:,3));
phi          = atan2d(GCOORD(:,2),GCOORD(:,1));
phi(phi < 0) = phi(phi < 0)+360;

GCOORD_SPH   = [theta phi r]; 

end % END OF FUNCTION cartesian2spherical