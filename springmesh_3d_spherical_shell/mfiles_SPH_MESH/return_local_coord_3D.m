function [lcoord] = return_local_coord_3D(gcoord,nodes,tetras,Pxyz)

% Part of 3D convection code M3TET_MG, which is a simplified version of the
% parallel 3D mantle convection code M3TET (developed by
% J.Hasenclever & J.Phipps Morgan, 2007-2010)
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH March 2011

ns = length(tetras); % Number of points
x  = reshape(gcoord(1,nodes(1:4,tetras)),4,ns);
y  = reshape(gcoord(2,nodes(1:4,tetras)),4,ns);
z  = reshape(gcoord(3,nodes(1:4,tetras)),4,ns);

xp = Pxyz(1,:);
yp = Pxyz(2,:);
zp = Pxyz(3,:);

lcoord = zeros(3,length(tetras));

lcoord(1,:) =   (  (x(2,:)-xp)    .*(y(3,:)-yp)    .*(z(4,:)-zp)     + (x(3,:)-xp)    .*(y(4,:)-yp)    .*(z(2,:)-zp)     + (x(4,:)-xp)    .*(y(2,:)-yp)    .*(z(3,:)-zp)     ...
                 - (z(2,:)-zp)    .*(y(3,:)-yp)    .*(x(4,:)-xp)     - (z(3,:)-zp)    .*(y(4,:)-yp)    .*(x(2,:)-xp)     - (z(4,:)-zp)    .*(y(2,:)-yp)    .*(x(3,:)-xp) )   ...
             ./ (  (x(2,:)-x(1,:)).*(y(3,:)-y(1,:)).*(z(4,:)-z(1,:)) + (x(3,:)-x(1,:)).*(y(4,:)-y(1,:)).*(z(2,:)-z(1,:)) + (x(4,:)-x(1,:)).*(y(2,:)-y(1,:)).*(z(3,:)-z(1,:))     ...
                 - (z(2,:)-z(1,:)).*(y(3,:)-y(1,:)).*(x(4,:)-x(1,:)) - (z(3,:)-z(1,:)).*(y(4,:)-y(1,:)).*(x(2,:)-x(1,:)) - (z(4,:)-z(1,:)).*(y(2,:)-y(1,:)).*(x(3,:)-x(1,:)) );

lcoord(2,:) =   (  (x(1,:)-xp)    .*(y(3,:)-yp)    .*(z(4,:)-zp)     + (x(3,:)-xp)    .*(y(4,:)-yp)    .*(z(1,:)-zp)     + (x(4,:)-xp)    .*(y(1,:)-yp)    .*(z(3,:)-zp)     ...
                 - (z(1,:)-zp)    .*(y(3,:)-yp)    .*(x(4,:)-xp)     - (z(3,:)-zp)    .*(y(4,:)-yp)    .*(x(1,:)-xp)     - (z(4,:)-zp)    .*(y(1,:)-yp)    .*(x(3,:)-xp) )   ...
             ./ (  (x(1,:)-x(2,:)).*(y(3,:)-y(2,:)).*(z(4,:)-z(2,:)) + (x(3,:)-x(2,:)).*(y(4,:)-y(2,:)).*(z(1,:)-z(2,:)) + (x(4,:)-x(2,:)).*(y(1,:)-y(2,:)).*(z(3,:)-z(2,:))     ...
                 - (z(1,:)-z(2,:)).*(y(3,:)-y(2,:)).*(x(4,:)-x(2,:)) - (z(3,:)-z(2,:)).*(y(4,:)-y(2,:)).*(x(1,:)-x(2,:)) - (z(4,:)-z(2,:)).*(y(1,:)-y(2,:)).*(x(3,:)-x(2,:)) );
             
lcoord(3,:) =   (  (x(1,:)-xp)    .*(y(2,:)-yp)    .*(z(4,:)-zp)     + (x(2,:)-xp)    .*(y(4,:)-yp)    .*(z(1,:)-zp)     + (x(4,:)-xp)    .*(y(1,:)-yp)    .*(z(2,:)-zp)     ...
                 - (z(1,:)-zp)    .*(y(2,:)-yp)    .*(x(4,:)-xp)     - (z(2,:)-zp)    .*(y(4,:)-yp)    .*(x(1,:)-xp)     - (z(4,:)-zp)    .*(y(1,:)-yp)    .*(x(2,:)-xp) )   ...
             ./ (  (x(1,:)-x(3,:)).*(y(2,:)-y(3,:)).*(z(4,:)-z(3,:)) + (x(2,:)-x(3,:)).*(y(4,:)-y(3,:)).*(z(1,:)-z(3,:)) + (x(4,:)-x(3,:)).*(y(1,:)-y(3,:)).*(z(2,:)-z(3,:))     ...
                 - (z(1,:)-z(3,:)).*(y(2,:)-y(3,:)).*(x(4,:)-x(3,:)) - (z(2,:)-z(3,:)).*(y(4,:)-y(3,:)).*(x(1,:)-x(3,:)) - (z(4,:)-z(3,:)).*(y(1,:)-y(3,:)).*(x(2,:)-x(3,:)) );

end