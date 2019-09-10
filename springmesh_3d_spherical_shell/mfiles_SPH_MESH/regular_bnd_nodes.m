function [pbnd] = regular_bnd_nodes(r,l0)
% Usage: [pbnd] = regular_bnd_nodes(r,l0)
%
% Purpose:
%   Create a regular distribution of boundary nodes in function of the 
%   distance between nodes (l0) on a sphere with radius r.
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

deg2rad   = pi/180;
ang_dist  = [31.68 15.79 7.96 3.98 1.99 0.99 0.50 0.25 0.13]';
d         = r*ang_dist*deg2rad;
[~,n]     = min(abs(d-l0));
% level   #triangles    #nodes	typ.angular dist
%    1           20         42	  31.68
%    2           80        162	  15.79
%    3          320        642	   7.96
%    4         1280       2562	   3.98
%    5         5120      10242	   1.99
%    6        20480      40962	   0.99
%    7        81920     163842     0.50
%    8       327680     655362     0.25
%    9      1310720    2621442     0.13
[pbnd,~] = sph_shell3(n,r);
pbnd     = pbnd';

end % END OF FUNCTION regular_bnd_nodes
