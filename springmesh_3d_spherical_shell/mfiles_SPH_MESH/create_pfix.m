function [pfix,pbnd] = create_pfix(pbnd,r)
% Usage: [pfix,pbnd] = create_pfix(pbnd,r)
%
% Purpose:
% Create 6 fixed nodes at the boundary (closest boundary nodes to theta =
% 0, 90, 180 and phi = 0, 90, 180 and 270 degrees) and then remove those
% nodes from pbnd. These pfix nodes give stability to the mesh when solving
% Hooke's law.
%
% Input:
%   pbnd  : [matrix]  : coordinates of boundary nodes
%   r     : [sccalar] : radius  
%
% Output:
%   pfix  : [matrix]  : coordinates of fixed nodes
%   pbnd  : [matrix]  : coordinates of boundary nodes after removing pfix 
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

pfix_ideal = [ 0  0  r ;...
               0  0 -r ;...
               r  0  0 ;...
               0  r  0 ;...
              -r  0  0 ;...
               0 -r  0];

index = zeros(1,6); % pre-allocate memory
for i = 1:6
    d = sqrt((pfix_ideal(i,1)-pbnd(:,1)).^2 + ...
             (pfix_ideal(i,2)-pbnd(:,2)).^2 + ...
             (pfix_ideal(i,3)-pbnd(:,3)).^2);
    [~,j]    = min(d);
    index(i) = j;
end
pfix          = pbnd(index,:); % create pfix
pbnd(index,:) = [];            % remove pfix values from pbnd

end % END OF SUBFUNCTION create_pfix