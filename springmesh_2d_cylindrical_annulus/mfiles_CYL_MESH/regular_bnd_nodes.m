function pbnd = regular_bnd_nodes(r,l0)
% Usage: pbnd = regular_bnd_nodes(r,l0)
%
% Purpose:
%   Create a regular distribution of boundary nodes in function of the 
%   distance between nodes (l0) on a circumference with radius r.
%
% Input:
%   r     : [scalar]        : circumference radius
%   l0    : [scalar]        : arc length (km) between boundary nodes
%
% Output:
%   pbnd  : [matrix]        : boundary nodes in Cartesian coordinates
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

deg2rad    = pi/180;
theta      = linspace(0,360,round(360/(l0/(r*deg2rad))))'; % Clokwise sense
theta(end) = [];                     % remove repeated angle, 0 = 360
pbnd       = zeros(length(theta),2); % pre-allocate memory
pbnd(:,1)  = r*sind(theta);          % x coordinate of points on the circumference of radius r
pbnd(:,2)  = r*cosd(theta);          % y coordinate of points on the circumference of radius r

end % END OF FUNCTION regular_bnd_nodes
