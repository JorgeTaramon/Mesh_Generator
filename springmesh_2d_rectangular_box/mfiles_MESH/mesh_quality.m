function q = mesh_quality(GCOORD,EL2NOD)
% Usage: q = mesh_quality(GCOORD,EL2NOD)
%
% Purpose:
%   Compute the quality factor for each triangle of the mesh
%
% Input: 
%   GCOORD : [matrix]        : npt-by-2 matrix with all the (x,y)
%                              coordinates 
%   EL2NOD : [matrix]        : nel-by-nnodel matrix, nnodel can be either 3
%                              or 6 (or 7) but we will only use the first 3
%                              cols (node 1,2,3) 
% Output:
%   q      : [column vector] : quality factor for each triangle
%
% Logic:
% "The plots of our 2-D meshes show that the algorithm produces triangles
% that are almost equilateral. This is a desirable property when solving
% PDEs with the Finite Element Method. Upper bounds on the errors depend
% only on the smallest angle in the mesh, and if all angles are close to 60
% degrees , good numerical results are achieved. The survey paper [24]
% discusses many measures of the "element quality". One commonly used
% quality measure is the ratio between the radius of the largest inscribed
% circle (times two) and the smallest circumscribed circle:
% q = 2 r_in/r_out = (b + c - a)(c + a - b)(a + b - c) / (a*b*c)   (2.29)
% where a, b, c are the side lengths. An equilateral triangle has q = 1,
% and a degenerate triangle (zero area) has q = 0. As a rule of thumb, if
% all triangles have q > 0.5 the results are good. "
%
%                                         --- cited from the DISTMESH paper
% Chao Jul 2008
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

bars   = [EL2NOD(:,[1,2]);EL2NOD(:,[2,3]);EL2NOD(:,[3,1])]; % we want to keep duplicated bars here
abc    = reshape(sqrt((GCOORD(bars(:,1),1)-GCOORD(bars(:,2),1)).^2 +...
                      (GCOORD(bars(:,1),2)-GCOORD(bars(:,2),2)).^2),size(EL2NOD,1),3);
q      = (abc(:,2) + abc(:,3) - abc(:,1)).*...
         (abc(:,3) + abc(:,1) - abc(:,2)).*...
         (abc(:,1) + abc(:,2) - abc(:,3))./...
         (abc(:,1).*abc(:,2).*abc(:,3));

%--------------- what the 'abc' and 'q' lines do --------------------------
% x_bar      = [GCOORD(bars(:,1),1) GCOORD(bars(:,2),1)];
% y_bar      = [GCOORD(bars(:,1),2) GCOORD(bars(:,2),2)];
% length_bar = sqrt((x_bar(:,1)-x_bar(:,2)).^2+(y_bar(:,1)-y_bar(:,2)).^2);
% nel        = size(EL2NOD,1);
% abc        = reshape(length_bar,nel,3);
% a          = abc(:,1);
% b          = abc(:,2);
% c          = abc(:,3);
% q          = (b + c - a).*(c + a - b).*(a + b - c) ./ (a.*b.*c);
%--------------------------------------------------------------------------

end % END OF FUNCTION mesh_quality