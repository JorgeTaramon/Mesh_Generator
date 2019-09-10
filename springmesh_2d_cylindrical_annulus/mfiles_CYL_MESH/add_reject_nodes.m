function [MESH] = add_reject_nodes(MESH,GUIDE_MESH,SETTINGS)
% Usage: [MESH] = add_reject_nodes(MESH,GUIDE_MESH,SETTINGS)
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
pbnd1      = MESH.pbnd1;
pbnd2      = MESH.pbnd2;
pint       = MESH.pint;
nfix       = size(MESH.pfix,1);
nbnd1      = size(MESH.pbnd1,1);
nbnd2      = size(MESH.pbnd2,1);
nint       = size(MESH.pint,1);
bars       = MESH.bars;
rel_change = MESH.rel_change;

if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    ifix1                   = find(abs(sqrt(sum(pfix.^2,2))-SETTINGS.r_int) < 1e-8); % indices for fixed nodes on boundary 1
    ibnd1                   = (nfix+1:nfix+nbnd1)';                                  % indices for bnd nodes on boundary 1
else
    ifix1                   = [];
    ibnd1                   = [];
end
ifix2                    = find(abs(sqrt(sum(pfix.^2,2))-SETTINGS.r_ext) < 1e-8); % indices for fixed nodes on boundary 2
ibnd2                    = (nfix+nbnd1+1:nfix+nbnd1+nbnd2)';                      % indices for bnd nodes on boundary 2
ifix_int                 = (1:nfix)';
ifix_int([ifix1; ifix2]) = [];                                                    % indices for fixed interior nodes

%==========================================================================
% COMPUTE NODES TO BE ADDED AND NODES TO BE REJECTED (BND1, BND2, INT)
%==========================================================================
bars_to_add_node        = bars(rel_change > 0.5,:); % stretched bars

i_bars_to_add_node_bnd1 = find(sum(ismember(bars_to_add_node,union(ifix1,ibnd1)),2) == 2); % bar indices to add a bnd1 node
i_bars_to_add_node_bnd2 = find(sum(ismember(bars_to_add_node,union(ifix2,ibnd2)),2) == 2); % bar indices to add a bnd2 node
bars_to_add_node_bnd1   = bars_to_add_node(i_bars_to_add_node_bnd1,:);          % bars to add a bnd1 node
bars_to_add_node_bnd2   = bars_to_add_node(i_bars_to_add_node_bnd2,:);          % bars to add a bnd2 node
bars_to_add_node_int    = bars_to_add_node;
bars_to_add_node_int([i_bars_to_add_node_bnd1;i_bars_to_add_node_bnd2],:) = []; % bars to add a int node

pbnd1_new = (GCOORD(bars_to_add_node_bnd1(:,1),:) + GCOORD(bars_to_add_node_bnd1(:,2),:))/2;
pbnd2_new = (GCOORD(bars_to_add_node_bnd2(:,1),:) + GCOORD(bars_to_add_node_bnd2(:,2),:))/2;
pint_new  = (GCOORD(bars_to_add_node_int(:,1),:)  + GCOORD(bars_to_add_node_int(:,2),:))/2;

%==========================================================================
% COMPUTE NODES TO BE REJECTED (BND1, BND2, INT)
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

i_bars_to_reject_node_bnd1 = find(sum(ismember(bars_to_reject_node_sorted,[ifix1;ifix_int;ibnd1]),2) == 2); % bar indices to reject a bnd1 node
i_bars_to_reject_node_bnd2 = find(sum(ismember(bars_to_reject_node_sorted,[ifix2;ifix_int;ibnd2]),2) == 2); % bar indices to reject a bnd2 node
bars_to_reject_node_bnd1   = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd1,:);   % bars to reject a bnd1 node
bars_to_reject_node_bnd2   = bars_to_reject_node_sorted(i_bars_to_reject_node_bnd2,:);   % bars to reject a bnd2 node
bars_to_reject_node_int    = bars_to_reject_node_sorted;
bars_to_reject_node_int([i_bars_to_reject_node_bnd1;i_bars_to_reject_node_bnd2],:) = []; % bars to reject a int node

bnd1_node_to_reject        = bars_to_reject_node_bnd1(:,2);
bnd2_node_to_reject        = bars_to_reject_node_bnd2(:,2);
int_node_to_reject         = bars_to_reject_node_int(:,2);

%===============================================================================
% REMOVE REPEATED NEW NODES AND PROJECT THE ADDED BND NODES ONTO THE BOUNDARIES 
%===============================================================================
% BOUNDARY 1
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    pbnd1_new = unique(pbnd1_new,'rows','stable');   % remove repeated pbnd1_new
    if ~isempty(pbnd1_new)
        % Project these new bnd1 nodes onto the circular boundary1
        pbnd1_new_POL        = cartesian2polar(pbnd1_new);                   % convert to polar coordinates
        pbnd1_new_POL(:,2)   = SETTINGS.r_int*ones(size(pbnd1_new_POL,1),1); % change the r-coordinate by r_int (projecting onto bnd 1)
        pbnd1_new            = polar2cartesian(pbnd1_new_POL);               % Cartesian coordinates of new bnd1 nodes after projecting on bnd1
        num_bnd1_nodes_added = size(pbnd1_new,1);                            % number of new bnd1 nodes
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
    % Project these new bnd2 nodes onto the circular boundary2
    pbnd2_new_POL        = cartesian2polar(pbnd2_new);                   % convert to polar coordinates
    pbnd2_new_POL(:,2)   = SETTINGS.r_ext*ones(size(pbnd2_new_POL,1),1); % change the r-coordinate by r_ext (projecting onto bnd 2)
    pbnd2_new            = polar2cartesian(pbnd2_new_POL);               % Cartesian coordinates of new bnd2 nodes after projecting on bnd2
    num_bnd2_nodes_added = size(pbnd2_new,1);                            % number of new bnd2 nodes
else
    num_bnd2_nodes_added = 0;
end
% INTERIOR
pint_new = unique(pint_new,'rows','stable'); % remove repeated pint_new
if ~isempty(pint_new)
    num_int_nodes_added = size(pint_new,1); % number of new int nodes
else
    num_int_nodes_added = 0;
end

%==========================================================================
% REMOVE REPEATED NODES TO REJECT
%==========================================================================
% BOUNDARY 1
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
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
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    pbnd1(bnd1_node_to_reject-nfix,:)       = []; % reject the 2nd node of the most compressed bar of the bad element (bnd1 nodes)
    pbnd1                                   = [pbnd1; pbnd1_new]; % add new bnd1 nodes
end
pbnd2(bnd2_node_to_reject-nfix-nbnd1,:)     = []; % reject the 2nd node of the most compressed bar of the bad element (bnd2 nodes)
pbnd2                                       = [pbnd2; pbnd2_new]; % add new bnd2 nodes
pint(int_node_to_reject-nfix-nbnd1-nbnd2,:) = []; % reject the 2nd node of the most compressed bar of the bad element (int nodes)
pint                                        = unique([pint; pint_new],'rows','stable'); % add new int nodes
GCOORD_new                                  = [pfix; pbnd1; pbnd2; pint];
EL2NOD_new                                  = delaunay(GCOORD_new);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_POL = cartesian2polar(GCOORD_new);
EL2NOD_new = EL2NOD_new(~(sum(ismember(EL2NOD_new,find(abs(GCOORD_POL(:,2)-SETTINGS.r_int) < 1e-8)),2)==3),:);
q_new      = mesh_quality(GCOORD_new,EL2NOD_new);

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd1                     = pbnd1;
MESH.pbnd2                     = pbnd2;
MESH.pint                      = pint;
MESH.GCOORD                    = GCOORD_new;
MESH.EL2NOD                    = EL2NOD_new;
MESH.q                         = q_new;
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,SETTINGS);
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
if SETTINGS.r_int > 0 && ~isempty(MESH.pbnd1)
    ADDREJECT.bnd1 = [num_bnd1_nodes_added num_bnd1_nodes_rejected];
end
ADDREJECT.bnd2 = [num_bnd2_nodes_added num_bnd2_nodes_rejected];
ADDREJECT.int  = [num_int_nodes_added  num_int_nodes_rejected];
display_add_reject_nodes(ADDREJECT)

end % END OF FUNCTION add_reject_nodes