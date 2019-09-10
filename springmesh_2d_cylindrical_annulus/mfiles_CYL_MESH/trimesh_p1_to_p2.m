function [GCOORD,EL2NOD,PointID] = trimesh_p1_to_p2(GCOORD_c,EL2NOD_c,PointID_c)
% Usage: [GCOORD,EL2NOD,PointID] = trimesh_p1_to_p2(GCOORD_c,EL2NOD_c,PointID_c)
% 
% Purpose: Returns a 6-node triangle connectivity matrix after adding edge
%          nodes to each element in a 3-node triangle mesh.
%
% Input:
%   GCOORD_c : [matrix] : coordinates of all nodes in 3-node triangle mesh
%   EL2NOD_c : [matrix] : finite element connectivity matrix (3 x nel)
%   PointID_c: [vector] : boundary indices for all nodes
%
% Output:
%   GCOORD : [matrix] : coordinates of all nodes in 6-node triangle mesh
%   EL2NOD : [matrix] : finite element connectivity matrix (6 x nel)
%   PointID: [vector] : boundary indices for all nodes
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
%

% Each linear 3-node element will get new edge nodes 4,5,6
% 
%  3
%  | \
%  6   5
%  |     \
%  1 - 4 - 2

% size of coarse mesh
nel_c  = size(EL2NOD_c,2);   % number of elements
nnod_c = size(GCOORD_c,2);   % number of nodes

% (1) Create quadratic order mesh of triangles
% =======================================================================
% Create a pointer edge to node defined by end-nodes
% (e.g. edge 1 has end-nodes [1 2], edge 2 has [2 3],...)
edge2node = reshape( EL2NOD_c([ 1  2
                               2  3
                               3  1 ]',:),2,[] )';

% Find the edges that are shared by neighboring elements and return a unique
% list of edges.
[edge2node,~,ib] = find_unique_edges(edge2node); % *SUBFUNCTION*
nedges  = size(edge2node,1); % number of unique edges (i.e. number of new edge 
                             % nodes that have to be generated)

EL2EGDE = reshape(ib,3,nel_c); % element to edge connectivity after double 
                               % edges have been removed

% Coordinates of the nodes defining each edge's end points
xBarEnds = reshape(GCOORD_c(1,edge2node'),2,[]);
yBarEnds = reshape(GCOORD_c(2,edge2node'),2,[]);

% Create new node at each edge mid-point
xBarMids = 0.5*sum(xBarEnds,1);
yBarMids = 0.5*sum(yBarEnds,1);

% storage for quadratic order (6-node) connectivity matrix
EL2NOD                  = zeros(6,nel_c,'uint32'); 
EL2NOD(1:3,:)           = EL2NOD_c;
EL2NOD(4:6,:)           = nnod_c+EL2EGDE;
% storage for quadratic order (6-node) node coordinates
nnod                    = nnod_c + nedges;
GCOORD                  = zeros(2,nnod);
GCOORD(:,1:nnod_c)      = GCOORD_c;
GCOORD(1,nnod_c+1:nnod) = xBarMids;
GCOORD(2,nnod_c+1:nnod) = yBarMids;

if nargin==3
    PointID           = zeros(1,nnod,'int32');
    PointID(1:nnod_c) = PointID_c;
    tmp               = PointID_c(edge2node);
    ind_0             = all(tmp,2)==0;
    ind_bnd           = zeros(nedges,1,'int32');
    ind_bnd(~ind_0)   = max(tmp(~ind_0,:),[],2);
    
    % Check for element edges that stretch from domain corner-to-corner
    % (they are domain edges, not corners)
    ind_c2c           = all(tmp,2)>100 & all(tmp,2)<200;
    if any(ind_c2c)
        error(' to be verified');
        % Node at bottom boundary
        ind           = all(ismember(tmp(ind_c2c,:),[101 102]),2);
        ind_bnd(ind2) = 201;
        % Node at right boundary
        ind           = all(ismember(tmp(ind_c2c,:),[102 103]),2);
        ind_bnd(ind2) = 202;
        % Node at top boundary
        ind           = all(ismember(tmp(ind_c2c,:),[103 104]),2);
        ind_bnd(ind2) = 203;
        % Node at left boundary
        ind           = all(ismember(tmp(ind_c2c,:),[101 104]),2);
        ind_bnd(ind2) = 204;
    end
    
    PointID(nnod_c+1:nnod) = ind_bnd;
end

end % END OF FUNCTION trimesh_p1_to_p2

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [edges_uniq,IIu,JJu] = find_unique_edges(edges)
    if size(edges,2)~=2
        error('expecting a nx2 matrix');
    end

    % find FIRST occurrence of unique rows in edges
    % use sort so that vertices defining the edges are sorted increasingly
    [tmp,iu,ju]   = unique(sort(edges,2),'rows','first');

    % create index IIu to extract edges_uniq from edges
    [IIu,tmp,JJu] = unique(iu(ju));

    % extract unique edges from input
    edges_uniq    = edges(IIu,:);
end % END OF SUBFUNCTION find_unique_edges
