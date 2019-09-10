function PointID = point_id_rectbox(GCOORD)
% Usage: PointID = point_id_rectbox(GCOORD)
% 
% Purpose: Generates indices for boundary nodes assuming rectangular a
%          domain.
%
% xmin        xmax          Coordinate system
% 104---203---103 zmax      0-------> x
%  |           |            |
%  |  domain   |            |
%  |           |            |
% 204         202           z
%  |           |
%  |           |
%  |           |
% 101---201---102  zmin
%
% Indices of domain boundaries
%
% Corners: 101 lower left
%          102 lower right
%          103 upper right
%          104 upper left
% Edges:   201 lower
%          202 right
%          203 upper
%          204 left
% Input:
%   GCOORD   : [matrix]    : Coordinates of all nodes in mesh
%
% Output:
%   PointID : [colvector] : node indices
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH : Dec 2012

tol  = 1e-8;
nnod = size(GCOORD,2);

xmax = max(GCOORD(1,:));
xmin = min(GCOORD(1,:));
zmax = max(GCOORD(2,:));
zmin = min(GCOORD(2,:));
Lx   = xmax-xmin;
Lz   = zmax-zmin;

% Right side nodes
nods_right = GCOORD(1,:)>xmax-tol*Lx;
% Left side nodes
nods_left  = GCOORD(1,:)<xmin+tol*Lx;
% Top side nodes
nods_top   = GCOORD(2,:)>zmax-tol*Lz;
% Bottom side nodes
nods_bot   = GCOORD(2,:)<zmin+tol*Lz;

PointID = zeros(1,nnod,'int32');
PointID(nods_bot)   = 201; % Bottom edge
PointID(nods_right) = 202; % Right edge
PointID(nods_top)   = 203; % Top edge
PointID(nods_left)  = 204; % Left edge
PointID(nods_left  & nods_bot) = 101; % Lower left corner
PointID(nods_right & nods_bot) = 102; % Lower right corner
PointID(nods_right & nods_top) = 103; % Upper right corner
PointID(nods_left  & nods_top) = 104; % Upper left corner
    
end % END OF FUNCTION point_id_rectbox