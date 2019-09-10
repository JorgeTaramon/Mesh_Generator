function coord = hexagonal_close_packing_spheres(xmin,xmax,ymin,ymax,zmin,zmax,r)
% Usage: coord = hexagonal_close_packing_spheres(xmin,xmax,ymin,ymax,zmin,zmax,r)
%
% Purpose:
%   Compute coordinates for the centres of spheres arranged in a simple hcp
%   lattice.
%
% Input:
%   xmin  : [scalar] : minimun x-coordinate
%   xmax  : [scalar] : maximum x-coordinate
%   ymin  : [scalar] : minimun y-coordinate
%   ymax  : [scalar] : maximum y-coordinate
%   zmin  : [scalar] : minimun z-coordinate
%   zmax  : [scalar] : maximum z-coordinate
%   r     : [scalar] : sphere radius
%
% Output:
%   coord : [matrix] : coorfinates for the spheres centres in a simple hcp
%                      lattice
%
% JMT Jan 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if nargin == 0
%     xmin = 0;
%     xmax = 4;
%     ymin = 0;
%     ymax = 3;
%     zmin = 0;
%     zmax = 2;
%     r    = 0.5;
    xmin = -6371;
    xmax =  6371;
    ymin = -6371;
    ymax =  6371;
    zmin = -6371;
    zmax =  6371;
    r    =  1000;
end

npt_x   = floor((xmax-xmin)/(2*r));
npt_y   = floor((ymax-ymin)/(sqrt(3)*r));
npt_z   = floor((zmax-zmin)/(2*sqrt(6)*r/3));
coord   = zeros(npt_x*npt_y*npt_z,3);
counter = 0;

for k = 0:npt_z-1
    for j = 0:npt_y-1
        for i = 0:npt_x-1
            coord(i+1+counter,1) = (2*i + rem(j+k,2))*r + xmin;
            coord(i+1+counter,2) = sqrt(3)*(j + rem(k,2)/3)*r + ymin;
            coord(i+1+counter,3) = (2/3)*sqrt(6)*k*r + zmin;
        end
        counter = counter + npt_x;
    end
end

% The function hexagonal_close_packing_spheres gives a compact hcp lattice
% but it starts to generate points from (xmin,ymin,zmin), so we do not have 
% a symmetric distribution of nodes from the origin (0,0,0). In order to do
% that, we have to center those coord.
mean_x      = mean(coord(:,1));
mean_y      = mean(coord(:,2));
mean_z      = mean(coord(:,3));
x_traslated = coord(:,1)-mean_x*ones(size(coord(:,1),1),1);
y_traslated = coord(:,2)-mean_y*ones(size(coord(:,2),1),1);
z_traslated = coord(:,3)-mean_z*ones(size(coord(:,3),1),1);
coord       = [x_traslated y_traslated z_traslated];

if nargin == 0 
    el2nod = delaunay(coord);
    q      = tetra_mesh_quality(coord,el2nod);
    figure(563);clf
    simpplot(coord,el2nod(q>=0.9,:))
    view(142.5,30)
    figure(564);clf
    simpplot(coord,el2nod(q>0.7 & q<0.9,:))
    view(142.5,30)
end

end % END OF SUBFUNCTION hexagonal_close_packing_spheres