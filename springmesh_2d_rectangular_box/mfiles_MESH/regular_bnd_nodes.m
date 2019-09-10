function pbnd = regular_bnd_nodes(segBC,l0,i_seg)
% Usage: pbnd = regular_bnd_nodes(segBC,l0,i_seg)
%
% Purpose:
%   Create a regular distribution of boundary nodes in a segment in
%   function of the distance between nodes (l0).
%
% Input:
%   segBC : [matrix] : info of boundary segments
%   l0    : [scalar] : length (km) between boundary nodes
%   i_seg : [scalar] : current segment
%
% Output:
%   pbnd  : [matrix] : boundary nodes in Cartesian coordinates
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

l_seg       = sqrt((segBC(i_seg,1) - segBC(i_seg,3)).^2 + (segBC(i_seg,2) - segBC(i_seg,4)).^2);
x           = linspace(segBC(i_seg,1),segBC(i_seg,3),round(l_seg/l0)+1);
z           = linspace(segBC(i_seg,2),segBC(i_seg,4),round(l_seg/l0)+1);
pbnd        = [x; z]';
pbnd(1,:)   = []; % remove pfix
pbnd(end,:) = []; % remove pfix

end % END OF FUNCTION regular_bnd_nodes
