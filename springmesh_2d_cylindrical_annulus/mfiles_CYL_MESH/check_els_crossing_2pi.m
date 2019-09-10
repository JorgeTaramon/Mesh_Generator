function [els_cross_2pi,els_no_cross_2pi] = check_els_crossing_2pi(GCOORD_POL,EL2NOD)
% Usage: [els_cross_2pi,els_no_cross_2pi] = check_els_crossing_2pi(GCOORD_POL,EL2NOD)
% 
% Purpose: Check which elements are crossing theta = 2pi and which ones are
%          not crossing theta = 2pi
%
% Input:
%   GCOORD_POL : [matrix] : polar coordinates for a 3/6/7/12-node triangle
%                           mesh 
%   EL2NOD     : [matrix] : finite element connectivity matrix (3 x nel)
%                           (6 x nel)(7 x nel)(12 x nel)
% Output:
%   els_cross_2pi    : [vector] : indices for those elements crossing
%                                 theta = 2pi 
%   els_no_cross_2pi : [vector] : indices for those elements not crossing
%                                 theta = 2pi
%
% Check if all nodes within the same element are separated a distance less
% than pi --> i.e. check if there is any element whose nodes can be in the
% 1st quadrant (0 < theta < 90) and in the 4th quadrant (270 < theta < 360)
% These nodes are close in Cartesian coordinates, however, they are ~2pi
% (360°) far from each other in polar coordinates.
% NOTE: theta (colatitude) is measured from +z axis in clockwise direction
% (range 0 to 2pi).
%
%                   Z
%                   |
%                   |
%                   |
%        1 *_- - - -|- - - - * 3
%             - _   |       /
%                 - |      /
%                   |- _  /
%                   |    * 2
%                   |
%                   |
%                   |
%   ----------------|----------------------> X
%                   |
%                   |
%                   |
%                   |
%                   |
%                   |
%
%
% JMT Oct 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

nel              = size(EL2NOD,2); % number of elements
VCOORD_th        = reshape(GCOORD_POL(1,EL2NOD(1:3,1:nel)),3,nel); % theta coordiante of the vertices
ang_dist         = zeros(3,nel); 
ang_dist(1,:)    = abs(VCOORD_th(1,:)-VCOORD_th(2,:)); % angular distance between vertex 1 and vertex 2
ang_dist(2,:)    = abs(VCOORD_th(2,:)-VCOORD_th(3,:)); % angular distance between vertex 2 and vertex 3
ang_dist(3,:)    = abs(VCOORD_th(3,:)-VCOORD_th(1,:)); % angular distance between vertex 3 and vertex 1
% Elements crossing theta = 2pi
els_cross_2pi = find(sum(ang_dist > pi,1)); % indices of those elements crossing 2pi
% Elements not crossing theta = 2pi
els_no_cross_2pi = (1:nel); % list of all elements
els_no_cross_2pi(ismember(els_no_cross_2pi,els_cross_2pi)) = []; % remove from the list those elements crossing 2pi

end % END OF FUNCTION check_els_crossing_2pi



