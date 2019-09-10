function coord = hexagonal_close_packing_circles(xmin,xmax,ymin,ymax,r)
% Usage: coord = hexagonal_close_packing_circles(xmin,xmax,ymin,ymax,r)
%
% Purpose:
%   Compute coordinates for the centres of circles arranged in a simple hcp
%   lattice. In this way we get regular triangles with edge length 2r.
%   It generates better triangles than using meshgrid built-in Matlab
%   function
%
% Input:
%   xmin  : [scalar] : minimun x-coordinate
%   xmax  : [scalar] : maximum x-coordinate
%   ymin  : [scalar] : minimun y-coordinate
%   ymax  : [scalar] : maximum y-coordinate
%   r     : [scalar] : circle radius
%
% Output:
%   coord : [matrix] : coorfinates for the centres of circles in a simple 
%                      hcp lattice
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

npt_x   = floor((xmax-xmin)/(2*r));
npt_y   = floor((ymax-ymin)/(sqrt(3)*r));
coord   = zeros(npt_x*npt_y,2);
counter = 0;

for j = 0:npt_y-1
    for i = 0:npt_x-1
        coord(i+1+counter,1) = (2*i + rem(j,2))*r + xmin;
        coord(i+1+counter,2) = sqrt(3)*j*r + ymin;
    end
    counter = counter + npt_x;
end

% % The function hexagonal_close_packing_circles gives a compact hcp lattice
% % but it starts to generate points from (xmin,ymin), so we do not have a 
% % symmetric distribution of nodes from the origin (0,0). In order to do
% % that, we have to center those coord.
% mean_x      = mean(coord(:,1));
% mean_y      = mean(coord(:,2));
% x_traslated = coord(:,1)-mean_x*ones(size(coord(:,1),1),1);
% y_traslated = coord(:,2)-mean_y*ones(size(coord(:,2),1),1);
% coord       = [x_traslated y_traslated];

end % END OF FUNCTION hexagonal_close_packing_circles