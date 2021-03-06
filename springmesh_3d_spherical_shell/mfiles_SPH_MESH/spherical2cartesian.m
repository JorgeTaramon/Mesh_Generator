function [GCOORD] = spherical2cartesian(GCOORD_SPH)
% Usage: [GCOORD] = spherical2cartesian(GCOORD_SPH)
%
% Purpose:
%   Transform spherical coordinates into Cartesian coordinates
%   
% Input:
%   GCOORD_SPH : [matrix] : Spherical coordinates (theta,phi,r) 
%                           (in degrees,degrees and km) 
%
% Output:
%   GCOORD     : [matrix] : Cartesian coordinates (x,y,z)(in km)
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

x = GCOORD_SPH(:,3).*sind(GCOORD_SPH(:,1)).*cosd(GCOORD_SPH(:,2));
y = GCOORD_SPH(:,3).*sind(GCOORD_SPH(:,1)).*sind(GCOORD_SPH(:,2));
z = GCOORD_SPH(:,3).*cosd(GCOORD_SPH(:,1));

GCOORD = [x y z];

end % END OF FUNCTION spherical2cartesian