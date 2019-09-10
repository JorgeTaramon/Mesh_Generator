function MESH = add_reject_nodes_v2(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: MESH = add_reject_nodes_v2(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   Add or reject nodes taking into account how is the relative bar-length
%   change which is given by:
%           Relative bar length change(i) = (L(i)-L0(i))/L0(i)
%   where L0 is the desired length and L is the actual length.
%   The criteria for either adding or rejecting nodes is:
%   - Adding a new node in those bars whose rel_change > 0.5, i.e., those
%     stretched bars more than a 50% regarding their desired length.
%   - Rejecting one of the nodes of the bars whose rel_change < -0.5, i.e.,
%     those compressed bars more than a 50% regarding their desired length.
%   This routine takes into account which kind of node (bnd1,bnd2 or int)
%   is being modified (either add or reject).
% 
% Input:
%   MESH       : [structure] : structure containing the mesh
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%   INTERFACE  : [structure] : structure containing interface settings
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   MESH       : [structure] : structure containing the mesh
%
% JMT Jan 2016
% JMT Jun 2016: cleaned up
% JMT Jan 2017: new add/reject version. This version is simpler and faster
%               than the old one (now called add_reject_nodes_old)
% JMT Aug 2017: Handles internal boundaries. Add/reject nodes taking into
%               account if they are on the internal boundaries.
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES AND COMPUTE NEEDED INFORMATION
%==========================================================================
r_int             = SETTINGS.r_int;
r_ext             = SETTINGS.r_ext;
GCOORD            = MESH.GCOORD;
pfix              = MESH.pfix;
pbnd1             = MESH.pbnd1;
pbnd2             = MESH.pbnd2;
pbnd_ref_face_bot = MESH.pbnd_ref_face_bot;
pbnd_ref_face_1   = MESH.pbnd_ref_face_1;
pbnd_ref_face_2   = MESH.pbnd_ref_face_2;
pbnd_ref_face_3   = MESH.pbnd_ref_face_3;
pbnd_ref_face_4   = MESH.pbnd_ref_face_4;
pint              = MESH.pint;
nfix              = size(pfix,1);
nbnd1             = size(pbnd1,1);
nbnd2             = size(pbnd2,1);
nbnd_ref_face_bot = size(pbnd_ref_face_bot,1);
nbnd_ref_face_1   = size(pbnd_ref_face_1,1);
nbnd_ref_face_2   = size(pbnd_ref_face_2,1);
nbnd_ref_face_3   = size(pbnd_ref_face_3,1);
nbnd_ref_face_4   = size(pbnd_ref_face_4,1);
nint              = size(pint,1);
bars              = MESH.bars;
rel_change        = MESH.rel_change;
deg2rad           = pi/180;
theta0            = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0              = GUIDE_MESH.phi0;     % longitude (degrees) of the point around which the refined and transition zones are defined
d_ref             = GUIDE_MESH.d_ref;    % refined zone depth (km)
w_ref_deg         = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees (North-South)
theta_ref_n       = theta0-w_ref_deg/2;  % colatitude of the northern boundary in the refined zone
theta_ref_s       = theta0+w_ref_deg/2;  % colatitude of the southern boundary in the refined zone
l_ref_deg         = GUIDE_MESH.l_ref/(deg2rad*r_ext);  % length of refined zone in degrees (East-West)
phi_ref_e         = phi0+l_ref_deg/2;    % longitude of the eastern boundary in the refined zone
phi_ref_w         = phi0-l_ref_deg/2;    % longitude of the western boundary in the refined zone
pfix_SPH          = cartesian2spherical(pfix);

if r_int > 0 && ~isempty(pbnd1)
    ifix1 = find(abs(pfix_SPH(:,3)-r_int) < 1e-8); % indices for fixed nodes on boundary 1
    ibnd1 = (nfix+1:nfix+nbnd1)';                  % indices for bnd nodes on boundary 1
else
    ifix1 = [];
    ibnd1 = [];
end
ifix2             = find(abs(pfix_SPH(:,3)-r_ext) < 1e-8);       % indices for fixed nodes on boundary 2
ibnd2             = (nfix+nbnd1+1:nfix+nbnd1+nbnd2)';            % indices for bnd nodes on boundary 2
ifix_ref_face_bot = find(abs(pfix_SPH(:,3)-(r_ext - d_ref)) < 1e-8); % indices for fixed nodes on bottom face of refined region
ibnd_ref_face_bot = (nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot)'; 
ifix_ref_face_1   = find(abs(pfix_SPH(:,1)-theta_ref_s) < 1e-8); % indices for fixed nodes on face 1 of refined region
ibnd_ref_face_1   = (nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)'; 
ifix_ref_face_2   = find(abs(pfix_SPH(:,2)-phi_ref_e) < 1e-8);   % indices for fixed nodes on face 2 of refined region
ibnd_ref_face_2   = (nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)'; 
ifix_ref_face_3   = find(abs(pfix_SPH(:,1)-theta_ref_n) < 1e-8); % indices for fixed nodes on face 3 of refined region
ibnd_ref_face_3   = (nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)'; 
ifix_ref_face_4   = find(abs(pfix_SPH(:,2)-phi_ref_w) < 1e-8);   % indices for fixed nodes on face 3 of refined region
ibnd_ref_face_4   = (nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)'; 
ifix_int          = (1:nfix)';
ifix_int([ifix1; ifix2;ifix_ref_face_bot;ifix_ref_face_1;ifix_ref_face_2;ifix_ref_face_3;ifix_ref_face_4]) = []; % indices for fixed interior nodes

%==========================================================================================
% COMPUTE NODES TO BE ADDED AND NODES TO BE REJECTED (BND1, BND2,REFINED REGION BNDS, INT)
%==========================================================================================
bars_to_add_node        = bars(rel_change > 0.5,:); % stretched bars

i_bars_to_add_node_bnd1             = find(sum(ismember(bars_to_add_node,union(ifix1,ibnd1)),2) == 2);                         % bar indices to add a bnd1 node
i_bars_to_add_node_bnd2             = find(sum(ismember(bars_to_add_node,union(ifix2,ibnd2)),2) == 2);                         % bar indices to add a bnd2 node
i_bars_to_add_node_bnd_ref_face_bot = find(sum(ismember(bars_to_add_node,union(ifix_ref_face_bot,ibnd_ref_face_bot)),2) == 2); % bar indices to add a bnd_ref_face_bot node
i_bars_to_add_node_bnd_ref_face_1   = find(sum(ismember(bars_to_add_node,union(ifix_ref_face_1,ibnd_ref_face_1)),2) == 2);     % bar indices to add a bnd_ref_face_1 node
i_bars_to_add_node_bnd_ref_face_2   = find(sum(ismember(bars_to_add_node,union(ifix_ref_face_2,ibnd_ref_face_2)),2) == 2);     % bar indices to add a bnd_ref_face_2 node
i_bars_to_add_node_bnd_ref_face_3   = find(sum(ismember(bars_to_add_node,union(ifix_ref_face_3,ibnd_ref_face_3)),2) == 2);     % bar indices to add a bnd_ref_face_3 node
i_bars_to_add_node_bnd_ref_face_4   = find(sum(ismember(bars_to_add_node,union(ifix_ref_face_4,ibnd_ref_face_4)),2) == 2);     % bar indices to add a bnd_ref_face_4 node

bars_to_add_node_bnd1               = bars_to_add_node(i_bars_to_add_node_bnd1,:);             % bars to add a bnd1 node
bars_to_add_node_bnd2               = bars_to_add_node(i_bars_to_add_node_bnd2,:);             % bars to add a bnd2 node
bars_to_add_node_bnd_ref_face_bot   = bars_to_add_node(i_bars_to_add_node_bnd_ref_face_bot,:); % bars to add a bnd_ref_face_bot node
bars_to_add_node_bnd_ref_face_1     = bars_to_add_node(i_bars_to_add_node_bnd_ref_face_1,:);   % bars to add a bnd_ref_face_1 node
bars_to_add_node_bnd_ref_face_2     = bars_to_add_node(i_bars_to_add_node_bnd_ref_face_2,:);   % bars to add a bnd_ref_face_2 node
bars_to_add_node_bnd_ref_face_3     = bars_to_add_node(i_bars_to_add_node_bnd_ref_face_3,:);   % bars to add a bnd_ref_face_3 node
bars_to_add_node_bnd_ref_face_4     = bars_to_add_node(i_bars_to_add_node_bnd_ref_face_4,:);   % bars to add a bnd_ref_face_4 node
bars_to_add_node_int                = bars_to_add_node;
bars_to_add_node_int([i_bars_to_add_node_bnd1;             ...
                      i_bars_to_add_node_bnd2;             ...
                      i_bars_to_add_node_bnd_ref_face_bot; ...
                      i_bars_to_add_node_bnd_ref_face_1;   ...
                      i_bars_to_add_node_bnd_ref_face_2;   ...
                      i_bars_to_add_node_bnd_ref_face_3;   ...
                      i_bars_to_add_node_bnd_ref_face_4],:) = []; % bars to add a int node

pbnd1_new             = (GCOORD(bars_to_add_node_bnd1(:,1),:)             + GCOORD(bars_to_add_node_bnd1(:,2),:))/2;
pbnd2_new             = (GCOORD(bars_to_add_node_bnd2(:,1),:)             + GCOORD(bars_to_add_node_bnd2(:,2),:))/2;
pbnd_ref_face_bot_new = (GCOORD(bars_to_add_node_bnd_ref_face_bot(:,1),:) + GCOORD(bars_to_add_node_bnd_ref_face_bot(:,2),:))/2;
pbnd_ref_face_1_new   = (GCOORD(bars_to_add_node_bnd_ref_face_1(:,1),:)   + GCOORD(bars_to_add_node_bnd_ref_face_1(:,2),:))/2;
pbnd_ref_face_2_new   = (GCOORD(bars_to_add_node_bnd_ref_face_2(:,1),:)   + GCOORD(bars_to_add_node_bnd_ref_face_2(:,2),:))/2;
pbnd_ref_face_3_new   = (GCOORD(bars_to_add_node_bnd_ref_face_3(:,1),:)   + GCOORD(bars_to_add_node_bnd_ref_face_3(:,2),:))/2;
pbnd_ref_face_4_new   = (GCOORD(bars_to_add_node_bnd_ref_face_4(:,1),:)   + GCOORD(bars_to_add_node_bnd_ref_face_4(:,2),:))/2;
pint_new              = (GCOORD(bars_to_add_node_int(:,1),:)              + GCOORD(bars_to_add_node_int(:,2),:))/2;

%==========================================================================
% COMPUTE NODES TO BE REJECTED (BND1, BND2, INT)
%==========================================================================
[rel_change_sorted,I]      = sort(rel_change);
I_nodes_to_be_rejected     = I(rel_change_sorted < -0.43);
bars_to_reject_node_sorted = bars(I_nodes_to_be_rejected,:); % bars sorted from the most compressed

% Find what bars have a node which has already been rejected, i.e., the node is in a bar situated over 
% those bars in the list bars_to_reject_node_sorted. This is done to avoid rejecting nodes that are linked 
[~,loc]                         = ismember(bars_to_reject_node_sorted(:,1),bars_to_reject_node_sorted(:,2));
iloc                            = find(loc>0);
bars_with_node_already_rejected = iloc(loc(loc>0) < iloc);
bars_to_reject_node_sorted(bars_with_node_already_rejected,:) = [];

i_bars_to_reject_node_bnd1             = find(sum(ismember(bars_to_reject_node_sorted,[ifix1;ifix_int;ibnd1]),2) == 2);                         % bar indices to reject a bnd1 node
i_bars_to_reject_node_bnd2             = find(sum(ismember(bars_to_reject_node_sorted,[ifix2;ifix_int;ibnd2]),2) == 2);                         % bar indices to reject a bnd2 node
i_bars_to_reject_node_bnd_ref_face_bot = find(sum(ismember(bars_to_reject_node_sorted,[ifix_ref_face_bot;ifix_int;ibnd_ref_face_bot]),2) == 2); % bar indices to reject a bnd_ref_face_bot node
i_bars_to_reject_node_bnd_ref_face_1   = find(sum(ismember(bars_to_reject_node_sorted,[ifix_ref_face_1;ifix_int;ibnd_ref_face_1]),2) == 2);     % bar indices to reject a bnd_ref_face_1 node
i_bars_to_reject_node_bnd_ref_face_2   = find(sum(ismember(bars_to_reject_node_sorted,[ifix_ref_face_2;ifix_int;ibnd_ref_face_2]),2) == 2);     % bar indices to reject a bnd_ref_face_2 node
i_bars_to_reject_node_bnd_ref_face_3   = find(sum(ismember(bars_to_reject_node_sorted,[ifix_ref_face_3;ifix_int;ibnd_ref_face_3]),2) == 2);     % bar indices to reject a bnd_ref_face_3 node
i_bars_to_reject_node_bnd_ref_face_4   = find(sum(ismember(bars_to_reject_node_sorted,[ifix_ref_face_4;ifix_int;ibnd_ref_face_4]),2) == 2);     % bar indices to reject a bnd_ref_face_4 node

bars_to_reject_node_bnd1               = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd1,:);             % bars to reject a bnd1 node
bars_to_reject_node_bnd2               = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd2,:);             % bars to reject a bnd2 node
bars_to_reject_node_bnd_ref_face_bot   = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd_ref_face_bot,:); % bars to reject a bnd_ref_face_bot node
bars_to_reject_node_bnd_ref_face_1     = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd_ref_face_1,:);   % bars to reject a bnd_ref_face_1 node
bars_to_reject_node_bnd_ref_face_2     = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd_ref_face_2,:);   % bars to reject a bnd_ref_face_2 node
bars_to_reject_node_bnd_ref_face_3     = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd_ref_face_3,:);   % bars to reject a bnd_ref_face_3 node
bars_to_reject_node_bnd_ref_face_4     = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd_ref_face_4,:);   % bars to reject a bnd_ref_face_4 node
bars_to_reject_node_int                = bars_to_reject_node_sorted;
bars_to_reject_node_int([i_bars_to_reject_node_bnd1;             ...
                         i_bars_to_reject_node_bnd2;             ...
                         i_bars_to_reject_node_bnd_ref_face_bot; ...
                         i_bars_to_reject_node_bnd_ref_face_1;   ...
                         i_bars_to_reject_node_bnd_ref_face_2;   ...
                         i_bars_to_reject_node_bnd_ref_face_3;   ...
                         i_bars_to_reject_node_bnd_ref_face_4],:) = []; % bars to add a int node

bnd1_node_to_reject             = bars_to_reject_node_bnd1(:,2);
bnd2_node_to_reject             = bars_to_reject_node_bnd2(:,2);
bnd_ref_face_bot_node_to_reject = bars_to_reject_node_bnd_ref_face_bot(:,2);
bnd_ref_face_1_node_to_reject   = bars_to_reject_node_bnd_ref_face_1(:,2);
bnd_ref_face_2_node_to_reject   = bars_to_reject_node_bnd_ref_face_2(:,2);
bnd_ref_face_3_node_to_reject   = bars_to_reject_node_bnd_ref_face_3(:,2);
bnd_ref_face_4_node_to_reject   = bars_to_reject_node_bnd_ref_face_4(:,2);
int_node_to_reject              = bars_to_reject_node_int(:,2);

%===============================================================================
% REMOVE REPEATED NEW NODES AND PROJECT THE ADDED BND NODES ONTO THE BOUNDARIES 
%===============================================================================
% BOUNDARY 1
if r_int > 0 && ~isempty(pbnd1)
    pbnd1_new = unique(pbnd1_new,'rows','stable'); % remove repeated pbnd1_new
    if ~isempty(pbnd1_new)
        % Project these new bnd1 nodes onto the spherical boundary1
        pbnd1_new_SPH        = cartesian2spherical(pbnd1_new);      % convert to spherical coordinates
        pbnd1_new_SPH(:,3)   = r_int*ones(size(pbnd1_new_SPH,1),1); % change the r-coordinate by r_int (projecting onto bnd 1)
        pbnd1_new            = spherical2cartesian(pbnd1_new_SPH);  % Cartesian coordinates of new bnd1 nodes after projecting on bnd1
        num_bnd1_nodes_added = size(pbnd1_new,1);                   % number of new bnd1 nodes
    else
        num_bnd1_nodes_added = 0;
    end
else
    pbnd1_new            = [];
    num_bnd1_nodes_added = 0;
end
% BOUNDARY 2
pbnd2_new = unique(pbnd2_new,'rows','stable'); % remove repeated pbnd2_new
if ~isempty(pbnd2_new)
    % Project these new bnd2 nodes onto the spherical boundary2
    pbnd2_new_SPH        = cartesian2spherical(pbnd2_new);      % convert to spherical coordinates
    pbnd2_new_SPH(:,3)   = r_ext*ones(size(pbnd2_new_SPH,1),1); % change the r-coordinate by r_ext (projecting onto bnd 2)
    pbnd2_new            = spherical2cartesian(pbnd2_new_SPH);  % Cartesian coordinates of new bnd2 nodes after projecting on bnd2
    num_bnd2_nodes_added = size(pbnd2_new,1);                   % number of new bnd2 nodes
else
    num_bnd2_nodes_added = 0;
end
% REFINE REGION BOTTOM FACE
pbnd_ref_face_bot_new = unique(pbnd_ref_face_bot_new,'rows','stable'); % remove repeated pbnd_ref_face_bot_new
if ~isempty(pbnd_ref_face_bot_new)
    % Project these new pbnd_ref_face_bot nodes onto the spherical boundary
    pbnd_ref_face_bot_new_SPH        = cartesian2spherical(pbnd_ref_face_bot_new);              % convert to spherical coordinates
    pbnd_ref_face_bot_new_SPH(:,3)   = (r_ext-d_ref)*ones(size(pbnd_ref_face_bot_new_SPH,1),1); % change the r-coordinate by r_ext-d_ref (projecting)
    pbnd_ref_face_bot_new            = spherical2cartesian(pbnd_ref_face_bot_new_SPH);          % Cartesian coordinates of new bnd_ref_face_bot nodes after projecting
    num_bnd_ref_face_bot_nodes_added = size(pbnd_ref_face_bot_new,1);                           % number of new bnd_ref_face_bot nodes
else
    num_bnd_ref_face_bot_nodes_added = 0;
end
% REFINE REGION FACE 1
pbnd_ref_face_1_new = unique(pbnd_ref_face_1_new,'rows','stable'); % remove repeated pbnd_ref_face_1_new
if ~isempty(pbnd_ref_face_1_new)
    num_bnd_ref_face_1_nodes_added = size(pbnd_ref_face_1_new,1);  % number of new bnd_ref_face_1 nodes
else
    num_bnd_ref_face_1_nodes_added = 0;
end
% REFINE REGION FACE 2
pbnd_ref_face_2_new = unique(pbnd_ref_face_2_new,'rows','stable'); % remove repeated pbnd_ref_face_2_new
if ~isempty(pbnd_ref_face_2_new)
    num_bnd_ref_face_2_nodes_added = size(pbnd_ref_face_2_new,1);  % number of new bnd_ref_face_2 nodes
else
    num_bnd_ref_face_2_nodes_added = 0;
end
% REFINE REGION FACE 3
pbnd_ref_face_3_new = unique(pbnd_ref_face_3_new,'rows','stable'); % remove repeated pbnd_ref_face_3_new
if ~isempty(pbnd_ref_face_3_new)
    num_bnd_ref_face_3_nodes_added = size(pbnd_ref_face_3_new,1);  % number of new bnd_ref_face_3 nodes
else
    num_bnd_ref_face_3_nodes_added = 0;
end
% REFINE REGION FACE 4
pbnd_ref_face_4_new = unique(pbnd_ref_face_4_new,'rows','stable'); % remove repeated pbnd_ref_face_4_new
if ~isempty(pbnd_ref_face_4_new)
    num_bnd_ref_face_4_nodes_added = size(pbnd_ref_face_4_new,1);  % number of new bnd_ref_face_4 nodes
else
    num_bnd_ref_face_4_nodes_added = 0;
end
% INTERIOR
pint_new = unique(pint_new,'rows','stable'); % remove repeated pint_new
if ~isempty(pint_new)
    num_int_nodes_added = size(pint_new,1);  % number of new int nodes
else
    num_int_nodes_added = 0;
end

%==========================================================================
% REMOVE REPEATED NODES TO REJECT
%==========================================================================
% BOUNDARY 1
if r_int > 0 && ~isempty(pbnd1)
    bnd1_node_to_reject = unique(bnd1_node_to_reject,'stable'); % remove repeated bnd1_node_to_reject
    if ~isempty(bnd1_node_to_reject)
        num_bnd1_nodes_rejected = size(bnd1_node_to_reject,1);  % number of bnd1 nodes to reject
    else
        num_bnd1_nodes_rejected = 0;
    end
else
    bnd1_node_to_reject     = [];
    num_bnd1_nodes_rejected = 0;
end
% BOUNDARY 2
bnd2_node_to_reject = unique(bnd2_node_to_reject,'stable'); % remove repeated bnd2_node_to_reject
if ~isempty(bnd2_node_to_reject)
    num_bnd2_nodes_rejected = size(bnd2_node_to_reject,1);  % number of bnd2 nodes to reject
else
    num_bnd2_nodes_rejected = 0;
end
% REFINE REGION BOTTOM FACE
bnd_ref_face_bot_node_to_reject = unique(bnd_ref_face_bot_node_to_reject,'stable'); % remove repeated pbnd_ref_face_bot_node_to_reject
if ~isempty(bnd_ref_face_bot_node_to_reject)
    num_bnd_ref_face_bot_nodes_rejected = size(bnd_ref_face_bot_node_to_reject,1);  % number of bnd_ref_face_bot nodes to reject
else
    num_bnd_ref_face_bot_nodes_rejected = 0;
end
% REFINE REGION FACE 1
bnd_ref_face_1_node_to_reject = unique(bnd_ref_face_1_node_to_reject,'stable'); % remove repeated pbnd_ref_face_1_node_to_reject
if ~isempty(bnd_ref_face_1_node_to_reject)
    num_bnd_ref_face_1_nodes_rejected = size(bnd_ref_face_1_node_to_reject,1);  % number of bnd_ref_face_1 nodes to reject
else
    num_bnd_ref_face_1_nodes_rejected = 0;
end
% REFINE REGION FACE 1
bnd_ref_face_2_node_to_reject = unique(bnd_ref_face_2_node_to_reject,'stable'); % remove repeated pbnd_ref_face_2_node_to_reject
if ~isempty(bnd_ref_face_2_node_to_reject)
    num_bnd_ref_face_2_nodes_rejected = size(bnd_ref_face_2_node_to_reject,1);  % number of bnd_ref_face_2 nodes to reject
else
    num_bnd_ref_face_2_nodes_rejected = 0;
end
% REFINE REGION FACE 3
bnd_ref_face_3_node_to_reject = unique(bnd_ref_face_3_node_to_reject,'stable'); % remove repeated pbnd_ref_face_3_node_to_reject
if ~isempty(bnd_ref_face_3_node_to_reject)
    num_bnd_ref_face_3_nodes_rejected = size(bnd_ref_face_3_node_to_reject,1);  % number of bnd_ref_face_3 nodes to reject
else
    num_bnd_ref_face_3_nodes_rejected = 0;
end
% REFINE REGION FACE 4
bnd_ref_face_4_node_to_reject = unique(bnd_ref_face_4_node_to_reject,'stable'); % remove repeated pbnd_ref_face_4_node_to_reject
if ~isempty(bnd_ref_face_4_node_to_reject)
    num_bnd_ref_face_4_nodes_rejected = size(bnd_ref_face_4_node_to_reject,1);  % number of bnd_ref_face_4 nodes to reject
else
    num_bnd_ref_face_4_nodes_rejected = 0;
end
% INTERIOR
int_node_to_reject  = unique(int_node_to_reject,'stable');  % remove repeated int_node_to_reject
if ~isempty(int_node_to_reject)
    num_int_nodes_rejected = size(int_node_to_reject,1);    % number of int nodes to reject
else
    num_int_nodes_rejected = 0;
end

%==========================================================================
% REJECT AND ADD THE NODES COMPUTED BEFORE
%==========================================================================
if r_int > 0 && ~isempty(pbnd1)
    pbnd1(bnd1_node_to_reject-nfix,:)       = []; % reject the 2nd node of the most compressed bar of the bad element (bnd1 nodes)
    pbnd1                                   = [pbnd1; pbnd1_new]; % add new bnd1 nodes
end
pbnd2(bnd2_node_to_reject-nfix-nbnd1,:)     = []; % reject the 2nd node of the most compressed bar of the bad element (bnd2 nodes)
pbnd2                                       = [pbnd2; pbnd2_new]; % add new bnd2 nodes
pbnd_ref_face_bot(bnd_ref_face_bot_node_to_reject-nfix-nbnd1-nbnd2,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd_ref_face_bot nodes)
pbnd_ref_face_bot                           = [pbnd_ref_face_bot; pbnd_ref_face_bot_new]; % add new bnd_ref_face_bot nodes
pbnd_ref_face_1(bnd_ref_face_1_node_to_reject-nfix-nbnd1-nbnd2-nbnd_ref_face_bot,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd_ref_face_1 nodes)
pbnd_ref_face_1                             = [pbnd_ref_face_1; pbnd_ref_face_1_new]; % add new bnd_ref_face_1 nodes
pbnd_ref_face_2(bnd_ref_face_2_node_to_reject-nfix-nbnd1-nbnd2-nbnd_ref_face_bot-nbnd_ref_face_1,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd_ref_face_2 nodes)
pbnd_ref_face_2                             = [pbnd_ref_face_2; pbnd_ref_face_2_new]; % add new bnd_ref_face_2 nodes
pbnd_ref_face_3(bnd_ref_face_3_node_to_reject-nfix-nbnd1-nbnd2-nbnd_ref_face_bot-nbnd_ref_face_1-nbnd_ref_face_2,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd_ref_face_3 nodes)
pbnd_ref_face_3                             = [pbnd_ref_face_3; pbnd_ref_face_3_new]; % add new bnd_ref_face_3 nodes
pbnd_ref_face_4(bnd_ref_face_4_node_to_reject-nfix-nbnd1-nbnd2-nbnd_ref_face_bot-nbnd_ref_face_1-nbnd_ref_face_2-nbnd_ref_face_3,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd_ref_face_4 nodes)
pbnd_ref_face_4                             = [pbnd_ref_face_4; pbnd_ref_face_4_new]; % add new bnd_ref_face_4 nodes
pint(int_node_to_reject-nfix-nbnd1-nbnd2-nbnd_ref_face_bot-nbnd_ref_face_1-nbnd_ref_face_2-nbnd_ref_face_3-nbnd_ref_face_4,:) = []; % reject the 2nd node of the most compressed bar of the bad element (int nodes)
pint                                        = unique([pint; pint_new],'rows','stable'); % add new int nodes
GCOORD_new                                  = [pfix;              ...
                                               pbnd1;             ...
                                               pbnd2;             ...
                                               pbnd_ref_face_bot; ...
                                               pbnd_ref_face_1;   ...
                                               pbnd_ref_face_2;   ...
                                               pbnd_ref_face_3;   ...
                                               pbnd_ref_face_4;   ...
                                               pint];
EL2NOD_new                                  = delaunay(GCOORD_new);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH  = cartesian2spherical(GCOORD_new);
EL2NOD_new  = EL2NOD_new(~(sum(ismember(EL2NOD_new,find(abs(GCOORD_SPH(:,3)-r_int) < 1e-8)),2)==4),:);
q_new       = tetra_mesh_quality(GCOORD_new,EL2NOD_new);

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd1                     = pbnd1;
MESH.pbnd2                     = pbnd2;
MESH.pbnd_ref_face_bot         = pbnd_ref_face_bot;
MESH.pbnd_ref_face_1           = pbnd_ref_face_1;
MESH.pbnd_ref_face_2           = pbnd_ref_face_2;
MESH.pbnd_ref_face_3           = pbnd_ref_face_3;
MESH.pbnd_ref_face_4           = pbnd_ref_face_4;
MESH.pint                      = pint;
MESH.GCOORD                    = GCOORD_new;
MESH.EL2NOD                    = EL2NOD_new;
MESH.q                         = q_new;
MESH.s                         = shape_measure(MESH.GCOORD,MESH.EL2NOD);
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
MESH.rel_change                = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
MESH.rel_change_abs            = abs(MESH.rel_change);        % absolute value of relative bar-length change
MESH.mean_misfit_bar_length    = sum(MESH.rel_change_abs)/size(MESH.L,1);
MESH.fraction_bnd1_changed     =    (num_bnd1_nodes_added + num_bnd1_nodes_rejected)/nbnd1;
MESH.net_change_bnd1           = abs(num_bnd1_nodes_added - num_bnd1_nodes_rejected)/nbnd1;
MESH.fraction_bnd2_changed     =    (num_bnd2_nodes_added + num_bnd2_nodes_rejected)/nbnd2;
MESH.net_change_bnd2           = abs(num_bnd2_nodes_added - num_bnd2_nodes_rejected)/nbnd2;
MESH.fraction_int_changed      =    (num_int_nodes_added  + num_int_nodes_rejected) /nint;
MESH.net_change_int            = abs(num_int_nodes_added  - num_int_nodes_rejected) /nint;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    figure(5)
    plot_net_change_nodes(MESH,SETTINGS)
end

%==========================================================================
% DISPLAY INFORMATION
%==========================================================================
if r_int > 0 && ~isempty(MESH.pbnd1)
    ADDREJECT.bnd1 = [num_bnd1_nodes_added num_bnd1_nodes_rejected];
end
ADDREJECT.bnd2             = [num_bnd2_nodes_added num_bnd2_nodes_rejected];
ADDREJECT.bnd_ref_face_bot = [num_bnd_ref_face_bot_nodes_added num_bnd_ref_face_bot_nodes_rejected];
ADDREJECT.bnd_ref_face_1   = [num_bnd_ref_face_1_nodes_added num_bnd_ref_face_1_nodes_rejected];
ADDREJECT.bnd_ref_face_2   = [num_bnd_ref_face_2_nodes_added num_bnd_ref_face_2_nodes_rejected];
ADDREJECT.bnd_ref_face_3   = [num_bnd_ref_face_3_nodes_added num_bnd_ref_face_3_nodes_rejected];
ADDREJECT.bnd_ref_face_4   = [num_bnd_ref_face_4_nodes_added num_bnd_ref_face_4_nodes_rejected];
ADDREJECT.int              = [num_int_nodes_added  num_int_nodes_rejected];
display_add_reject_nodes(ADDREJECT)

end % END OF FUNCTION add_reject_nodes_v2