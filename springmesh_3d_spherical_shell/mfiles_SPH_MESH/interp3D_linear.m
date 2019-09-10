function V_i = interp3D_linear(nodes,els,lcoord,V)
%
% Purpose:
% Linearly interpolate variables at given local coordinates in elements
%
% Input
%   nodes  :: element connectivity
%   els    :: element in which each interpolation point is located
%   lcoord :: local coordinates of each interpolation point in "els"
%   V      :: variable that is about to be interpolated
% Output
%   V_i    :: interpolated values of V
%
% Part of 3D convection code M3TET_MG, which is a simplified
% version of the 3D mantle convection code M3TET (developed by
% J.Hasenclever & J.Phipps Morgan, 2007-2010)
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH April 2011

verbose = 0;
if verbose
    disp   ('************************************************')
    fprintf('interpol_lin: range lcoord %4.2e %4.2e %4.2e %4.2e\n',...
             min(lcoord(1,:)), max(lcoord(1,:)), min(lcoord(2,:)), max(lcoord(2,:)));
    disp   ('************************************************')
end

% Calculate linear shape function values at local coordinates
N      = [lcoord;1-sum(lcoord)];
N(N>1) = 1; % To avoid strange interpolation values due to round off errors
N(N<0) = 0;

% Interpolate to points using shape functions
V_i    = sum(N.*V(nodes(1:4,els)));
V_i    = V_i(:);

end % END OF FUNCTION interp3D_linear