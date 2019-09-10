function L0 = bar_L0_interface(GCOORD,h0,INTERFACE)
% Usage: L0 = bar_L0_interface(GCOORD,h0,INTERFACE)
%
% Purpose:
%   Compute L0 values for each node (GCOORD) in a narrow region given by a
%   spherical interface (based on Anderson et al., 2005).
%
% Input:
%   GCOORD     : [matrix]        : points in which L0 will be computed
%   h0         : [scalar]        : desired spring length (km) for a regular
%                                  mesh 
%   INTERFACE  : [structure]     : structure containing interface settings
%
% Output:
%   L0         : [column vector] : desired length for the springs
%
% JMT Jan 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

barmid_sph = cartesian2spherical(GCOORD);
L1         = h0*ones(size(GCOORD,1),1);
L2         = L1/INTERFACE.ref_factor + INTERFACE.s*(abs(barmid_sph(:,3)-repmat(INTERFACE.r,size(GCOORD,1),1)));

L0 = min(L1,L2);

end % END OF FUNCTION bar_L0_interface