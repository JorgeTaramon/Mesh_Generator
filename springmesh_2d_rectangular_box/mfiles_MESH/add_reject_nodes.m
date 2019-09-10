function MESH = add_reject_nodes(MESH,GUIDE_MESH,SETTINGS)
% Usage: MESH = add_reject_nodes(MESH,GUIDE_MESH,SETTINGS)
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
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   MESH       : [structure] : structure containing the mesh
%
% JMT May 2016
% JMT Jun 2017: new add/reject version. This version is simpler and faster
%               than the old one (now called add_reject_nodes_old)
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES AND COMPUTE NEEDED INFORMATION
%==========================================================================
GCOORD     = MESH.GCOORD;
pfix       = MESH.pfix;
pbnd       = MESH.pbnd;
pint       = MESH.pint;
nfix       = size(MESH.pfix,1);
nbnd       = size(MESH.pbnd,1);
nint       = size(MESH.pint,1);
bars       = MESH.bars;
rel_change = MESH.rel_change;
ifix       = (1:nfix)';
ibnd       = (nfix+1:nfix+nbnd)';
pseg       = MESH.pseg;

%==========================================================================
% COMPUTE NODES TO BE ADDED AND NODES TO BE REJECTED (BND, INT)
%==========================================================================
bars_to_add_node       = bars(rel_change > 0.5,:); % stretched bars

i_bars_to_add_node_bnd = find(sum(ismember(bars_to_add_node,union(ifix,ibnd)),2) == 2); % bar indices to add a bnd node
bars_to_add_node_bnd   = bars_to_add_node(i_bars_to_add_node_bnd,:);          % bars to add a bnd node
bars_to_add_node_int   = bars_to_add_node;
bars_to_add_node_int(i_bars_to_add_node_bnd,:) = []; % bars to add a int node

if strcmp(SETTINGS.refinement,'guide_mesh')
    nfix_mesh_vertex = size(MESH.pfix_mesh_vertex,1);
    ifix_mesh_vertex = (1:nfix_mesh_vertex)';
    ifix_ref_vertex  = nfix_mesh_vertex + (1:size(MESH.pfix_ref_vertex,1))';
    ifix_ref_vertex(MESH.pfix_ref_vertex(:,2) == GUIDE_MESH.z0 - GUIDE_MESH.d_ref) = []; % remove from the list vertices on the bottom of refined region
    ifix_mesh_vertex = [ifix_mesh_vertex; ifix_ref_vertex];
    
    i_bars_to_add_node_bnd = find(sum(ismember(bars_to_add_node,union(ifix_mesh_vertex,ibnd)),2) == 2); % bar indices to add a bnd node
    bars_to_add_node_bnd   = bars_to_add_node(i_bars_to_add_node_bnd,:);                         % bars to add a bnd node
    bars_to_add_node_int   = bars_to_add_node;
    bars_to_add_node_int(i_bars_to_add_node_bnd,:) = []; % bars to add a int node
end

pbnd_new = (GCOORD(bars_to_add_node_bnd(:,1),:) + GCOORD(bars_to_add_node_bnd(:,2),:))/2;
pint_new = (GCOORD(bars_to_add_node_int(:,1),:) + GCOORD(bars_to_add_node_int(:,2),:))/2;

pseg_new = pseg(bars_to_add_node_bnd(:,2)-nfix); % segment_ID for added bnd nodes

%==========================================================================
% REMOVE REPEATED NEW NODES
%==========================================================================
% BOUNDARY
pbnd_new = unique(pbnd_new,'rows','stable'); % remove repeated pbnd_new
if ~isempty(pbnd_new)
    num_bnd_nodes_added = size(pbnd_new,1);  % number of new bnd nodes
else
    num_bnd_nodes_added = 0;
end

% INTERIOR
pint_new = unique(pint_new,'rows','stable'); % remove repeated pint_new
if ~isempty(pint_new)
    num_int_nodes_added = size(pint_new,1);  % number of new int nodes
else
    num_int_nodes_added = 0;
end

%==========================================================================
% COMPUTE NODES TO BE REJECTED (BND, INT)
%==========================================================================
[rel_change_sorted,I]      = sort(rel_change);
I_nodes_to_be_rejected     = I(rel_change_sorted < -0.5);
bars_to_reject_node_sorted = bars(I_nodes_to_be_rejected,:); % bars sorted from the most compressed

% Find what bars have a node which has already been rejected, i.e., the node is in a bar situated over 
% those bars in the list bars_to_reject_node_sorted. This is done to avoid rejecting nodes that are linked 
[~,loc]                         = ismember(bars_to_reject_node_sorted(:,1),bars_to_reject_node_sorted(:,2));
iloc                            = find(loc>0);
bars_with_node_already_rejected = iloc(loc(loc>0) < iloc);
bars_to_reject_node_sorted(bars_with_node_already_rejected,:) = [];

i_bars_to_reject_node_bnd  = find(sum(ismember(bars_to_reject_node_sorted,[ifix;ibnd]),2) == 2); % bar indices to reject a bnd node
bars_to_reject_node_bnd    = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd,:);   % bars to reject a bnd node
bars_to_reject_node_int    = bars_to_reject_node_sorted;
bars_to_reject_node_int(i_bars_to_reject_node_bnd,:) = []; % bars to reject a int node

bnd_node_to_reject         = bars_to_reject_node_bnd(:,2);
int_node_to_reject         = bars_to_reject_node_int(:,2);

%==========================================================================
% REMOVE REPEATED NODES TO REJECT
%==========================================================================
% BOUNDARY
bnd_node_to_reject = unique(bnd_node_to_reject,'stable'); % remove repeated bnd_node_to_reject
if ~isempty(bnd_node_to_reject)
    num_bnd_nodes_rejected = size(bnd_node_to_reject,1);  % number of bnd nodes to reject
else
    num_bnd_nodes_rejected = 0;
end

% INTERIOR
int_node_to_reject = unique(int_node_to_reject,'stable'); % remove repeated int_node_to_reject
if ~isempty(int_node_to_reject)
    num_int_nodes_rejected = size(int_node_to_reject,1);  % number of int nodes to reject
else
    num_int_nodes_rejected = 0;
end

%==========================================================================
% REJECT AND ADD THE NODES COMPUTED BEFORE
%==========================================================================
pbnd(bnd_node_to_reject-nfix,:)      = []; % reject the 2nd node of the most compressed bar of the bad element (bnd nodes)
pseg(bnd_node_to_reject-nfix)        = []; % remove the segment_ID for the rejected bnd nodes
pbnd                                 = [pbnd; pbnd_new]; % add new bnd nodes
pseg                                 = [pseg, pseg_new]; % add segment_ID for added bnd nodes
pint(int_node_to_reject-nfix-nbnd,:) = []; % reject the 2nd node of the most compressed bar of the bad element (int nodes)
pint                                 = unique([pint; pint_new],'rows','stable'); % add new int nodes
GCOORD_new                           = [pfix; pbnd; pint];
EL2NOD_new                           = delaunay(GCOORD_new);
q_new                                = mesh_quality(GCOORD_new,EL2NOD_new);

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd                      = pbnd;
MESH.pint                      = pint;
MESH.GCOORD                    = GCOORD_new;
MESH.EL2NOD                    = EL2NOD_new;
MESH.pseg                      = pseg;
MESH.q                         = q_new;
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,SETTINGS);
MESH.rel_change                = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
MESH.rel_change_abs            = abs(MESH.rel_change);        % absolute value of relative bar-length change
MESH.mean_misfit_bar_length    = sum(MESH.rel_change_abs)/size(MESH.L,1);
MESH.fraction_bnd_changed      =    (num_bnd_nodes_added + num_bnd_nodes_rejected)/nbnd;
MESH.net_change_bnd            = abs(num_bnd_nodes_added - num_bnd_nodes_rejected)/nbnd;
MESH.fraction_int_changed      =    (num_int_nodes_added + num_int_nodes_rejected) /nint;
MESH.net_change_int            = abs(num_int_nodes_added - num_int_nodes_rejected) /nint;

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
ADDREJECT.bnd = [num_bnd_nodes_added num_bnd_nodes_rejected];
ADDREJECT.int = [num_int_nodes_added num_int_nodes_rejected];
display_add_reject_nodes(ADDREJECT)

end % END OF FUNCTION add_reject_nodes