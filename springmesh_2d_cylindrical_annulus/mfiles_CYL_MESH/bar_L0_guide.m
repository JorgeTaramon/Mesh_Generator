function L0 = bar_L0_guide(GCOORD,GUIDE_MESH)
% Usage: L0 = bar_L0_guide(GCOORD,GUIDE_MESH)
%
% Purpose:
%   Compute L0 values for each node (GCOORD) interpolating over the guide
%   mesh in polar coordinates (much more accurate than in Cartesian
%   coordinates because of the straight edges)
%
% Input:
%   GCOORD     : [matrix]        : points in which L0 will be computed
%   GUIDE_MESH : [structure]     : structure containing guide mesh
%
% Output:
%   L0         : [column vector] : interpolated values of desired length
%                                  for the springs
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

L0_guide     = GUIDE_MESH.L0_guide;
GCOORD_GUIDE = GUIDE_MESH.GCOORD_GUIDE';
EL2NOD_GUIDE = uint32(GUIDE_MESH.EL2NOD_GUIDE');
GCOORD_POL   = cartesian2polar(GCOORD);
GCOORD_POL   = GCOORD_POL';
r_int        = GUIDE_MESH.r_int;
r_ext        = GUIDE_MESH.r_ext;
r_tol        = 1e-12 * (r_ext - r_int); % Tolerance to make sure that bnd nodes are exactly on the boundary 
                                        % (the radius of bnd nodes can shift slightly (~1e-12) from the actual
                                        % boundary due to coordinate transformation)
GCOORD_POL(2,GCOORD_POL(2,:) <= r_int + r_tol) = r_int;
GCOORD_POL(2,GCOORD_POL(2,:) >= r_ext - r_tol) = r_ext;

% For every node (x,y) of GCOORD, find which element (EL2NOD_GUIDE) is in
els          = tsearch2(GCOORD_GUIDE,EL2NOD_GUIDE,GCOORD_POL); 
% Find local coordinates (r,s) which correspond to (x,y) in each element
local_coord  = local_coords_2d(GCOORD_GUIDE,EL2NOD_GUIDE,els,GCOORD_POL);
% Find out L0 in GCOORD_POL interpolating the values of L0_guide in local 
% coordinates through shape functions
L0           = interp2d_tri367(EL2NOD_GUIDE,els,local_coord,L0_guide);

end % END OF FUNCTION bar_L0_guide