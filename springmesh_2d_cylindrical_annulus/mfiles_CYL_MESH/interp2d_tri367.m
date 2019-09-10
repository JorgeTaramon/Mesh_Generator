function V_i = interp2d_tri367(EL2NOD,els,lc,V)
% Usage: V_i = interp2d_tri367(EL2NOD,els,lc,V)
%
% Purpose: Interpolates variables in a triangular finite element mesh using
%          linear or quadratic interpolation functions.
%
% Input
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : elements in which will be interpolated
%   lc     : [matrix]    : local coordinates in each element
%   V      : [colvector] : variable field to be interpolated
%
% Output
%   V_i    : [colvector] : interpolated values
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JH Dec 2012 : function will interpolate linearly (3-node triangles) and
%               quadratically (6- and 7-node triangles) depending on size
%               of EL2NOD
% JH Jan 2013 : : handles NaN in lc

ind    = ~isnan(lc(1,:)) & ~isnan(lc(2,:));
els    = els(ind);

% Local coordinates
r      = lc(1,ind);
s      = lc(2,ind);
nnodel = size(EL2NOD,1);
N      = sf_tri367(nnodel,r,s); % *SUBFUNCTION*

% Interpolate temperature to points using shape functions
V_i      = nan(length(ind),1);
V_i(ind) = sum(N'.*V(EL2NOD(:,els)));

end % END OF FUNCTION interp2d_tri367

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function N = sf_tri367(nnodel,r,s)
% Find shape functions and their derivatives at given points on the
% master element for a 7 node triangle

r   = r(:);
s   = s(:);
t   = 1 - r - s;
npt = length(r);
N   = zeros(npt,nnodel);
switch nnodel
    case 3
        % 3-node triangle
        %
        %        3
        %        | \
        % s-axis |   \
        %        |     \
        %        1 - - - 2     
        %          r axis -->
        %
        N(:,1) = 1 - r - s; % N1 at coordinate (r,s)
        N(:,2) = r;         % N2 at coordinate (r,s)
        N(:,3) = s;         % N3 at coordinate (r,s)
        
    case 6
        % 6-node triangle (node numbering is important)
        %
        %        3
        %        | \
        % s-axis 6   5
        %        |     \
        %        1 - 4 - 2
        %          r-axis
        %
        N(:,1) = t.*(2.*t-1) ; % N1 at coordinate (r,s)
        N(:,2) = r.*(2.*r-1) ; % N2 at coordinate (r,s)
        N(:,3) = s.*(2.*s-1) ; % etc
        N(:,4) = 4.*t.*r     ;
        N(:,5) = 4.*r.*s     ;
        N(:,6) = 4.*t.*s     ;
        
    case 7
        % 7-node triangle
        %
        %        3
        %        | \
        % s-axis 6   5
        %        | 7   \
        %        1 - 4 - 2
        %          r-axis
        %
        % There are 2 different versions for a 7-node triangle:
        % VERSION (1) has 7 shape functions that sum up to 1 everywhere in
        % the triangle.
        % VERSION (2) has 6 shape functions identical to those of a standard
        % 6-node triangle plus a 'bubble' centroid function that is 1 in the
        % center and zero at the other 6 nodes. Variables associated with
        % the bubble functions are DEVIATIONS of velocities from the
        % quadratic solution of a 6-node triangle.
        % IMPORTANT: Make sure you use the same 7-node version here and in
        %            function "shp_deriv_tri367"!!
        
        % VERSION (1)
        N(:,1) = t.*(2.*t-1)+ 3.*t.*r.*s;
        N(:,2) = r.*(2.*r-1)+ 3.*t.*r.*s;
        N(:,3) = s.*(2.*s-1)+ 3.*t.*r.*s;
        N(:,4) = 4.*r.*s - 12.*t.*r.*s  ;
        N(:,5) = 4.*t.*s - 12.*t.*r.*s  ;
        N(:,6) = 4.*t.*r - 12.*t.*r.*s  ;
        N(:,7) = 27.*t.*r.*s            ;
        
%         % VERSION (2)
%         % FIRST 6 functions IDENTICAL TO 6-node triangles
%         N(:,1) = t.*(2.*t-1) ; % N1 at coordinate (r,s)
%         N(:,2) = r.*(2.*r-1) ; % N2 at coordinate (r,s)
%         N(:,3) = s.*(2.*s-1) ; % etc
%         N(:,4) = 4.*t.*r     ;
%         N(:,5) = 4.*r.*s     ;
%         N(:,6) = 4.*t.*s     ;
%         % BUBBLE with 1 at center (1/3,1/3,1/3), 0 at nodes 1-6
%         N(:,7) = 27.*r.*s.*t ; 
end

end % END OF FUNCTION sf_tri367