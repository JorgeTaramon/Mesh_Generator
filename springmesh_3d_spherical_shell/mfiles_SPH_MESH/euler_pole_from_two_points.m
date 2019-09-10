function [EP_lat,EP_lon,EP_angle] = euler_pole_from_two_points(P1_SPH,P2_SPH)
% Usage: [EP_lat,EP_lon,EP_angle] = euler_pole_from_two_points(P1_SPH,P2_SPH)
%
% Purpose:
%   Compute the Euler Pole (lat,lon) and angle rotated between 2 points
%
% Input:
%   P1_SPH   : [vector] : spherical coordinates for point 1
%   P2_SPH   : [vector] : spherical coordinates for point 2
%
% Output:
%   EP_lat   : [scalar] : Euler Pole latitude (degrees)
%   EP_lon   : [scalar] : Euler Pole longitude (degrees)
%   EP_angle : [scalar] : rotated angle (degrees)
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

P1       = spherical2cartesian(P1_SPH);
P2       = spherical2cartesian(P2_SPH);
C        = cross(P1,P2);
C_SPH    = cartesian2spherical(C);
EP_lat   = 90 - C_SPH(1,1);                       % Euler Pole latitude
EP_lon   = C_SPH(1,2);                            % Euler Pole longitude
EP_angle = acosd(dot(P1,P2)/(norm(P1)*norm(P2))); % Rotation angle

end % END OF FUNCTION euler_pole_from_two_points