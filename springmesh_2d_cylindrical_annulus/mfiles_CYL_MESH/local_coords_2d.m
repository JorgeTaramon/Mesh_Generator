function LCOORD = local_coords_2d(GCOORD,EL2NOD,els,gX_PT)
% Usage: LCOORD = local_coords_2d(GCOORD,EL2NOD,els,gX_PT)
%
% Purpose: Returns local coordinates of points "gX_PT" in elements "els".
%
% Input:
%   GCOORD : [matrix]    : coordinates of all nodes in mesh
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : element in which each point is located
%   gX_PT  : [matrix]    : coordinates of points to be located
% Output:
%   LCOORD : [matrix]    : local coordinates of points in each element
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JH Jan 2013 : handles NaN in els

ind = ~isnan(els) & els>0;
els = els(ind);
ns  = length(els);
x   = reshape(GCOORD(1,EL2NOD(1:3,els)),3,ns);
y   = reshape(GCOORD(2,EL2NOD(1:3,els)),3,ns);

xp  = gX_PT(1,ind);
yp  = gX_PT(2,ind);

LCOORD = nan(2,length(ind));
LCOORD(1,ind) = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:).*yp-xp.*y(3,:))./ ...
                 (-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));
LCOORD(2,ind) = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp.*y(1,:)-xp.*y(2,:))./ ...
               (-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));

end % END OF FUNCTION local_coords_2d