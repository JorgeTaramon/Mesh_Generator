function s = shape_measure(GCOORD,EL2NOD)
% Usage: s = shape_measure(GCOORD,EL2NOD)
%
% Purpose:
%   Compute the shape measure for each tetrahedron of the mesh. The 
%   shape_measure is given by (Anderson et al., 2005):
%
%                       12(9V^2)^(1/3)    
%               s = -----------------------
%                    sum_i=1_to_6((h_i)^2)
% 
% Input: 
%   GCOORD : [matrix]        : npt-by-3 matrix with all the (x,y,z)
%                              coordinates 
%   EL2NOD : [matrix]        : nel-by-nnodel matrix, nnodel can be either 4
%                              or 10 but we will only use the first 4 cols
%                              (node 1,2,3,4) 
% Output:
%   s      : [column vector] : shape measure for each tetrahedron
%
% JMT Jan 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

V      = calc_tetra_volume(GCOORD',EL2NOD');
V      = -1.*V; % The connectivity for computing the tetrahedron volume which follows Hughes' notation
                % and the connectivity given by Delaunay are different. One is [1 2 3 4] and the other
                % one [1 2 4 3]. For that reason we have to multiply by -1.
if any(V(V<0) < -1e-6)
    error('no good')
end

bars   = [EL2NOD(:,[1,2]); ...
          EL2NOD(:,[1,3]); ...
          EL2NOD(:,[1,4]); ...
          EL2NOD(:,[2,3]); ...
          EL2NOD(:,[2,4]); ...
          EL2NOD(:,[3,4])];
barvec = GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:); % list of bar vectors
L      = reshape(sqrt(sum(barvec.^2,2)),size(EL2NOD,1),6); % bar lengths for each element

s      = (12*(9*V.^2).^(1/3))./sum(L.^2,2);

end % END OF FUNCTION shape_measure