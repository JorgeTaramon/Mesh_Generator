function [GCOORD_curved,GCOORD_POL,EL2NOD,PointID] = ...
    trimesh_p1_to_p2_cyl(GCOORD_c,GCOORD_POL_c,EL2NOD_c,PointID_c)
% Usage: [GCOORD_curved,GCOORD_POL,EL2NOD,PointID] = ...
%   trimesh_p1_to_p2_cyl(GCOORD_c,GCOORD_POL_c,EL2NOD_c,PointID_c)
% 
% Purpose: Returns a 6-node triangle connectivity matrix after adding edge
%          nodes to each element in a 3-node triangle mesh. This routine
%          returns 6-node triangles with curved egdes since it is splitting
%          the edges in polar coordinates. Remember that straight egde in
%          polar coordinates (theta,r) means curved edge in cartesian
%          coordinates (x,z).
%
% Input:
%   GCOORD_c      : [matrix] : Cartesian coordinates of all nodes in 3-node
%                              triangle mesh
%   GCOORD_POL_c  : [matrix] : polar coordinates of all nodes in 3-node
%                              triangle mesh
%   EL2NOD_c      : [matrix] : finite element connectivity matrix (3 x nel)
%   PointID_c     : [vector] : boundary indices for all nodes
%
% Output:
%   GCOORD_curved        : [matrix] : Cartesian coordinates of all nodes 
%                                     new mesh (2 x ?) 
%   GCOORD_POL           : [matrix] : polar coordinates of all nodes in 
%                                     6-node triangle mesh
%   EL2NOD               : [matrix] : finite element connectivity matrix 
%                                     (6 x nel)
%   PointID              : [vector] : boundary indices for all nodes
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JMT Mar 2016 : Compute curved edges for a polar mesh (theta,r)
% JMT Oct 2016 : handles elements crossing theta = 2pi with a 180° rotation
%                instead of using extended matrices (old version)
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%======================================================================================================
% COMPUTE ELEMENTS CROSSING THETA = 2PI
%======================================================================================================
[els_cross_2pi,~] = check_els_crossing_2pi(GCOORD_POL_c,EL2NOD_c);

%======================================================================================================
% COMPUTE MID-SIDE NODES IN POLAR COORDINATES FOR ELEMENTS 
% CROSSING THETA = 2PI AFTER A 180° ROTATION AROUND Y AXIS 
%======================================================================================================
[GCOORD_POL_Mids_els_cross_2pi,bar2node_els_cross_2pi] = mid_side_nodes_pol_rot_frame_els_cross_2pi...
    (GCOORD_c,EL2NOD_c,els_cross_2pi);

%======================================================================================================
% COMPUTE MID-SIDE NODES IN POLAR COORDINATES FOR ALL ELEMENTS
%======================================================================================================
[GCOORD_curved,GCOORD_POL,EL2NOD,PointID] = mid_side_nodes_pol_frame_all_els...
    (GCOORD_POL_c,EL2NOD_c,GCOORD_POL_Mids_els_cross_2pi,bar2node_els_cross_2pi,PointID_c);

end % END OF FUNCTION trimesh_p1_to_p2_cyl

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [edges_uniq,IIu,JJu] = find_unique_bars(edges)
    if size(edges,2)~=2
        error('expecting a nx2 matrix');
    end

    % find FIRST occurrence of unique rows in edges
    % use sort so that vertices defining the edges are sorted increasingly
    [~,iu,ju]   = unique(sort(edges,2),'rows','first');

    % create index IIu to extract edges_uniq from edges
    [IIu,~,JJu] = unique(iu(ju));

    % extract unique edges from input
    edges_uniq    = edges(IIu,:);
end % END OF SUBFUNCTION find_unique_edges

% #########################################################################

function [GCOORD_POL_Mids_els_cross_2pi,bar2node_els_cross_2pi] = ...
    mid_side_nodes_pol_rot_frame_els_cross_2pi(GCOORD_c,EL2NOD_c,els_cross_2pi)
% Usage: [GCOORD_POL_Mids_els_cross_2pi,bar2node_els_in_cone] = ...
%   mid_side_nodes_pol_rot_frame_els_cross_2pi(GCOORD_c,EL2NOD_c,els_cross_2pi)
%
% Purpose: 
%   Compute mid-side nodes for elements crossing theta = 2pi in a rotated
%   polar frame (180° around Y axis) following the next steps:
%       - Make a 180° counterclockwise rotation around the Y axis for the
%         elements crossing theta = 2pi.
%       - Convert to polar coordinates.
%       - Compute mid-side nodes in straight edges of the rotated polar
%         frame.
%       - Convert back to Cartesian coordinates.
%       - Undo the rotation (make a 180° clockwise rotation around the Y
%         axis).
%       - Convert again to polar coordinates: GCOORD_POL_Mids_els_cross_2pi
%
% Input:
%   GCOORD_c      : [matrix] : Cartesian coordinates of all nodes (3 nodel)
%                              (2 x nnod_c)
%   EL2NOD_c      : [matrix] : connectivity matrix (3 x nel)
%   els_cross_2pi : [vector] : indices for those elements crossing
%                              theta = 2pi 
%
% Output:
%   GCOORD_POL_Mids_els_cross_2pi : [matrix] : midside nodes in polar  
%                                              coordinates for elements 
%                                              crossing theta = 2pi
%   bar2node_els_cross_2pi        : [matrix] : bar list for elements 
%                                              crossing theta = 2pi
%
% JMT Oct 2016

% Make a 180° counterclockwise rotation around the Z axis for elements crossing theta = 2pi
RR_Z_180_CCW                      = [-1  0 ; ...
                                      0 -1 ];
GCOORD_c_rot                      = RR_Z_180_CCW * GCOORD_c;

% Convert to polar coordinates
GCOORD_POL_c_rot                  = cartesian2polar_rad(GCOORD_c_rot);

% Compute mid-side nodes in straight edges of the polar system
EL2NOD_c_els_cross_2pi            = EL2NOD_c(:,els_cross_2pi);
bar2node_els_cross_2pi            = reshape( EL2NOD_c_els_cross_2pi([ 1  2 ... % 1st edge of parent
                                                                      2  3 ... % 3nd edge of parent
                                                                      3  1 ... % 3rd edge of parent
                                                                           ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars
[bar2node_els_cross_2pi,~,~]      = find_unique_bars(bar2node_els_cross_2pi);

% Coordinates of the nodes defining each bar's end points
thetaBarEnds_els_cross_2pi        = reshape(GCOORD_POL_c_rot(1,bar2node_els_cross_2pi'),2,[]);
rBarEnds_els_cross_2pi            = reshape(GCOORD_POL_c_rot(2,bar2node_els_cross_2pi'),2,[]);

% Create new node at each bar mid-point
thetaBarMids_els_cross_2pi        = 0.5*sum(thetaBarEnds_els_cross_2pi,1);
rBarMids_els_cross_2pi            = 0.5*sum(rBarEnds_els_cross_2pi,1);
GCOORD_POL_Mids_els_cross_2pi_rot = [thetaBarMids_els_cross_2pi; ...
                                     rBarMids_els_cross_2pi];
 
% Convert to Cartesian coordinates
GCOORD_Mids_els_cross_2pi_rot     = polar2cartesian_rad(GCOORD_POL_Mids_els_cross_2pi_rot);

% Undo the rotation (rotate 180° around Y axis clockwise)
RR_Y_180_CW                       = [-1  0 ; ...
                                      0 -1 ];
GCOORD_Mids_els_cross_2pi         = RR_Y_180_CW * GCOORD_Mids_els_cross_2pi_rot;

% Convert to polar coodinates
GCOORD_POL_Mids_els_cross_2pi     = cartesian2polar_rad(GCOORD_Mids_els_cross_2pi);

end % END OF SUBFUNCTION mid_side_nodes_pol_rot_frame_els_cross_2pi

% #########################################################################

function [GCOORD_curved,GCOORD_POL,EL2NOD,PointID] = mid_side_nodes_pol_frame_all_els... 
    (GCOORD_POL_c,EL2NOD_c,GCOORD_POL_Mids_els_cross_2pi,bar2node_els_cross_2pi,PointID_c)
% Usage: [GCOORD_curved,GCOORD_POL,EL2NOD,PointID] = mid_side_nodes_pol_frame_all_els... 
%   (GCOORD_POL_c,EL2NOD_c,GCOORD_POL_Mids_els_cross_2pi,bar2node_els_cross_2pi,PointID_c)
%
% Purpose: 
%   Compute mid-side nodes for all elements in the original polar frame 
%   using the mid-side nodes of elements crossing theta = 2pi computed in
%   mid_side_nodes_pol_rot_frame_els_cross_2pi subfunction.
%   The steps are:  
%       - Compute bars for all elements (bar2node).
%       - Compute mid-side nodes in straight edges of the original polar
%         frame.
%       - Find bar positions of "bar2node_els_cross_2pi" in "bar2node".
%       - Substitute GCOORD_POL_Mids_els_cross_2pi in the mid-side nodes 
%         computed in straight edges of the original polar frame using the 
%         bar positions of "bar2node_els_cross_2pi" in "bar2node".
%
% Input:
%   GCOORD_POL_c                  : [matrix] : polar coordinates of all 
%                                              nodes (3 nodel) (2 x nnod_c)
%   EL2NOD_c                      : [matrix] : connectivity matrix
%                                              (3 x nel)
%   GCOORD_POL_Mids_els_cross_2pi : [matrix] : midside nodes in polar  
%                                              coordinates for elements 
%                                              crossing theta = 2pi
%   bar2node_els_cross_2pi        : [matrix] : bar list for elements 
%                                              crossing theta = 2pi
%   PointID_c                     : [vector] : boundary indices for all 
%                                              nodes
%
% Output:
%   GCOORD_curved : [matrix] : Cartesian coordinates of all nodes new mesh
%                              (2 x ?) 
%   GCOORD_POL    : [matrix] : polar coordinates of all nodes new mesh
%                              (2 x ?) 
%   EL2NOD        : [matrix] : connectivity matrix (6 x nel)
%   PointID       : [vector] : boundary indices for all nodes
%
% JMT Oct 2016

% Size of coarse mesh
nel_c           = size(EL2NOD_c,2);     % number of elements
nnod_c          = size(GCOORD_POL_c,2); % number of nodes

% Create a pointer bars defined by their end-nodes (e.g. bar 1 has end-nodes [1 5], bar 2 has [5 2],...)
bar2node        = reshape( EL2NOD_c([ 1  2 ... % 1st edge of parent
                                      2  3 ... % 3nd edge of parent
                                      3  1 ... % 3rd edge of parent
                                           ]',:),2,[] )';

% Find the bars that are shared by neighboring elements and return a unique list of bars
[bar2node,~,ib] = find_unique_bars(bar2node); % *SUBFUNCTION*

nbars           = size(bar2node,1);    % number of unique bars
EL2BAR          = reshape(ib,3,nel_c); % element to bar connectivity after doubles have been merged

% Coordinates of the nodes defining each bar's end points
thetaBarEnds    = reshape(GCOORD_POL_c(1,bar2node'),2,[]);
rBarEnds        = reshape(GCOORD_POL_c(2,bar2node'),2,[]);

% Create new node at each bar mid-point
thetaBarMids    = 0.5*sum(thetaBarEnds,1);
rBarMids        = 0.5*sum(rBarEnds,1);

% Check if any new node in phiBarMids has theta > 2pi, if so subtract 2pi
thetaBarMids(thetaBarMids > 2*pi) = thetaBarMids(thetaBarMids > 2*pi) - 2*pi;

% Compute positions of bar2node_els_cross_2pi in the bar2node list
[~,LOCB]  = ismember(bar2node_els_cross_2pi,bar2node,'rows');
[~,LOCB2] = ismember(fliplr(bar2node_els_cross_2pi),bar2node,'rows'); % this is just in case some bars were flipped
                                                                      % when using find_unique_bars subfunction
LOCB(LOCB == 0) = LOCB2(LOCB2 ~= 0); clear LOCB2;
% Substitute mid-side nodes of elements crossing theta = 2pi
thetaBarMids(LOCB) = GCOORD_POL_Mids_els_cross_2pi(1,:);
rBarMids(LOCB)     = GCOORD_POL_Mids_els_cross_2pi(2,:);

% Storage for quadratic order (10-node) connectivity matrix
EL2NOD        = zeros(6,nel_c,'uint32'); 
EL2NOD(1:3,:) = EL2NOD_c;
EL2NOD(4:6,:) = nnod_c+EL2BAR;

% Storage for quadratic order (10-node) node coordinates
nnod                        = nnod_c + nbars;
GCOORD_POL                  = zeros(2,nnod);
GCOORD_POL(:,1:nnod_c)      = GCOORD_POL_c;
GCOORD_POL(1,nnod_c+1:nnod) = thetaBarMids;
GCOORD_POL(2,nnod_c+1:nnod) = rBarMids;

% Convert to Cartesian coordinates
GCOORD_curved = polar2cartesian_rad(GCOORD_POL);

if nargin==5
    PointID           = zeros(1,nnod,'int32');
    PointID(1:nnod_c) = PointID_c;
    tmp               = PointID_c(bar2node);
    ind_0             = all(tmp,2)==0;
    ind_bnd           = zeros(nbars,1,'int32');
    ind_bnd(~ind_0)   = max(tmp(~ind_0,:),[],2);
    PointID(nnod_c+1:nnod) = ind_bnd;
end

end % END OF SUBFUNCTION mid_side_nodes_pol_frame_all_els