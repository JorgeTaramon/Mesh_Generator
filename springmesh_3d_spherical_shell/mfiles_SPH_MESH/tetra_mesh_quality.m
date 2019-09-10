function q = tetra_mesh_quality(GCOORD,EL2NOD)
% Usage: q = tetra_mesh_quality(GCOORD,EL2NOD)
%
% Purpose:
%   Compute the quality factor for each tetrahedron of the mesh. The 
%   quality factor is given by:
%
%                          Radius inscribed sphere   
%               q = 3 * -----------------------------
%                        Radius circumscribed sphere 
%
%   An equilateral tetrahredron has the max quality factor, q = 1
%   A degenerated tetrahedron (coplanar) has the min quality factor, q = 0
%
%   If OABC forms a generalized tetrahedron with a vertex O as the origin
%   and vectors a, b, and c, represent the positions of the vertices A, B
%   and C with respect to O, then the radius of the insphere is given by:
%
%           r = 6V /(|bxc| + |cxa| +|axb| + |(bxc)+(cxa)+(axb)|)
%
%   and the radius of the circumsphere is given by:
%
%           R = | a^2(bxc) + b^2(cxa) + c^2(axb) | / (12V)
%
%   where
%
%           6V = |a . (b x c)|
%
% Input: 
%   GCOORD : [matrix]        : npt-by-3 matrix with all the (x,y,z)
%                              coordinates 
%   EL2NOD : [matrix]        : nel-by-nnodel matrix, nnodel can be either 4
%                              or 10 but we will only use the first 4 cols
%                              (node 1,2,3,4) 
% Output:
%   q      : [column vector] : quality factor for each tetrahedron
%
% Chao Mar 2009
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

a     = GCOORD(EL2NOD(:,2),:) - GCOORD(EL2NOD(:,1),:);
b     = GCOORD(EL2NOD(:,3),:) - GCOORD(EL2NOD(:,1),:);
c     = GCOORD(EL2NOD(:,4),:) - GCOORD(EL2NOD(:,1),:);
axb   = cross(a,b,2);
bxc   = cross(b,c,2);
cxa   = cross(c,a,2);
tempR = repmat(dot(a,a,2),1,3).*bxc + ...
        repmat(dot(b,b,2),1,3).*cxa + ...
        repmat(dot(c,c,2),1,3).*axb;    
SixV  = dot(a,bxc,2);
r     = SixV ./ ( sqrt(sum(bxc.*bxc,2)) + ...
                  sqrt(sum(cxa.*cxa,2)) + ...
                  sqrt(sum(axb.*axb,2)) + ...
                  sqrt(sum((bxc+cxa+axb).*(bxc+cxa+axb),2))); % radius of inscribed sphere
R     = sqrt((sum(tempR.*tempR,2)))./SixV/2; % radius of circumscribed sphere
q     = 3*r./R; % quality factor

end % END OF FUNCTION tetra_mesh_quality