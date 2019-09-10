function MESH = fix_slivers(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: MESH = fix_slivers(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   Find slivers (q < q_sliver) and fix them by adding or rejecting some
%   nodes. A sliver is an almost coplanar tetrahedron which quality factor
%   is very close to zero. This routine is purely geometric (it does not
%   take into account the desired or actual length of the springs), and
%   replaces a sliver by the best triangle that can be formed using some of
%   its nodes.
%   For each sliver, compute the possible combinations of 3 nodes that give
%   the best triangle (best quality factor).
%   For example, being the sliver with nodes 1,2,3 and 4: 
%                                           
%                           4                            
%                           .                            
%                        .  |     .                      
%                     .     |          .                 
%                  .         |             .            
%               1. - - - - - - - - - - - - - -. 3       
%                  .          |            .             
%                     .       |         .                 
%                        .     |     .                   
%                           .  |  .                      
%                              .                        
%                              2                     
%                                              
%   there are 8 or 9 possible triangles associated to:
%   - 1 vertex and 2 midpoints (CASE 1):
%       Triangle formed by nodes 1, midpoint in 2-3 and modipoint in 3-4
%       Triangle formed by nodes 2, midpoint in 1-4 and modipoint in 3-4
%       Triangle formed by nodes 3, midpoint in 1-2 and modipoint in 1-4
%       Triangle formed by nodes 4, midpoint in 1-2 and modipoint in 2-3
%   - 3 vertices (1 tetrahedon face) (CASE 2):
%       Triangle formed by nodes 1, 2 and 3
%       Triangle formed by nodes 1, 2 and 4
%       Triangle formed by nodes 1, 3 and 4
%       Triangle formed by nodes 2, 3 and 4
%   
%   If there is also a shared edge, e.g. 1-2, there is another possible
%   triangle associated to: 
%   - 2 vertices and 1 midpoint (of the shared bar)(CASE 3):
%       Triangle formed by midpoint in 1-2, node 3 and node 4
%
%   Then, quality factor for these 8 or 9 triangles is compared, selecting
%   the one with the best quality factor. If e.g. the best triangle is the
%   one formed by nodes 4, midpoint in 1-2 and modipoint in 2-3 the
%   routine:
%       * Add midpoint between 1 and 2, and midpoint between 2 and 3
%       * Remove nodes 1, 2 and 3
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
% JMT Nov 2015
% JMT Jul 2016: cleaned up. Now boundary nodes in slivers can be added or
%               removed
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% DEBUG MODE
%==========================================================================
%
if nargin == 0
    SETTINGS.q_sliver                   = 0.1;
    SETTINGS.r_int                      = [];
    SETTINGS.r_ext                      = [];
    SETTINGS.h0                         = 1;
    SETTINGS.mean_misfit_bar_length_tol = 0.14;
    SETTINGS.guide_mesh                 = 'no';
    SETTINGS.show_figs                  = 1;
    SETTINGS.save_figs                  = 0;
    GUIDE_MESH                          = struct([]);
    INTERFACE                           = struct([]);
    MESH                                = choose_case(SETTINGS);
end

%==========================================================================
% LOAD VARIABLES
%==========================================================================
q_sliver             = SETTINGS.q_sliver;
q                    = MESH.q;
EL2NOD               = MESH.EL2NOD;
GCOORD               = MESH.GCOORD;
pfix                 = MESH.pfix;
pbnd1                = MESH.pbnd1;
pbnd2                = MESH.pbnd2;
pint                 = MESH.pint;
nfix                 = size(pfix,1);
nbnd1                = size(pbnd1,1);
nbnd2                = size(pbnd2,1);
pbnd1_new            = []; % empty pbnd1_new to enter in the loop
pbnd2_new            = []; % empty pbnd2_new to enter in the loop
pint_new             = []; % empty pint_new to enter in the loop
bnd1_nodes_to_remove = []; % empty bnd1_node_to_remove to enter in the loop
bnd2_nodes_to_remove = []; % empty bnd2_node_to_remove to enter in the loop
int_nodes_to_remove  = []; % empty int_node_to_remove to enter in the loop

%==========================================================================
% FIND SLIVERS AND FIX THEM
%==========================================================================
SLIVER2NOD        = EL2NOD(q < q_sliver,:);         % elements with q < q_sliver
[~,iq]            = sort(q(q < q_sliver),'ascend'); % sort quality factor of the slivers from the worst one to the best one
SLIVER2NOD_sorted = SLIVER2NOD(iq,:);               % sort the slivers from the worst one
bars_all_slivers  = sort(transpose(reshape(transpose( ...
                         SLIVER2NOD_sorted(:,[1 2     ...
                                              1 3     ...
                                              1 4     ...
                                              2 3     ...
                                              2 4     ...
                                              3 4])),2,[])),2); % sliver bars (also the repeated ones)

% NOTE: A 'for' loop is used in order to deal with silvers individually
% since some of them could share some nodes and bars, and this means an
% issue when we deal with all slivers at the same time. This implementation
% ('for' loop) is not very time consuming since there are rarely too many
% slivers. 
for i = 1:size(SLIVER2NOD_sorted,1)
    sorted_nodes  = sort(SLIVER2NOD_sorted(i,:),2); % sort nodes for this sliver from smallest one to the biggest one 
                                                    % (this is done in order to find the sliver bars easily)
    I_fix         = find(sorted_nodes <= nfix); % local index of fix nodes for this sliver
    num_bnd_nodes = sum( sorted_nodes > nfix & ...
                         sorted_nodes <= nfix + nbnd1 + nbnd2); % number of bnd nodes for this sliver
    if SETTINGS.r_int > 0 && ~isempty(pbnd1)
        I_bnd1    = find(sorted_nodes >  nfix                 & ...
                         sorted_nodes <= nfix + nbnd1);         % local index of bnd1 nodes for this sliver
    else
        I_bnd1    = [];
    end
    I_bnd2        = find(sorted_nodes >  nfix + nbnd1         &...
                         sorted_nodes <= nfix + nbnd1 + nbnd2); % local index of bnd2 nodes for this sliver
    I_int         = find(sorted_nodes >  nfix + nbnd1 + nbnd2); % local index of int nodes for this sliver
    bars          = transpose(reshape(transpose(  ...
                              sorted_nodes(1,[1 2 ...
                                              1 3 ...
                                              1 4 ...
                                              2 3 ...
                                              2 4 ...
                                              3 4])),2,[])); % bars for this sliver sorted in two columns 
                                                             % (each colum is the node at each end of the bar)
    barmid        = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoints of the bars for this sliver
    % Check if any node of this sliver has already been removed because of it was shared with some other previous sliver in the loop: 
    %   - In case that none of the nodes have been removed, proceed normally
    %   - In case that any node has been removed, do nothing more and go to the next sliver in the loop (since this sliver has already been fixed)
    sorted_nodes_no_fixed_nodes        = sorted_nodes;
    sorted_nodes_no_fixed_nodes(I_fix) = []; % remove fixed nodes from the list (since fixed nodes cannot be removed) 
    check_fixed_sliver = sum(ismember(sorted_nodes_no_fixed_nodes', [bnd1_nodes_to_remove; bnd2_nodes_to_remove; int_nodes_to_remove])) == 0;
    if check_fixed_sliver % check_fixed_sliver = 1 --> no node has been removed
                          % check_fixed_sliver = 0 --> at least one node has been removed
        %==============================================================================================================
        % FIND OUT IF THE SLIVER IS SHARING SOME EDGES (shared_bars) WITH OTHER SLIVERS AND HOW MANY (num_shared_bars)
        %==============================================================================================================
        [~,ind]         = ismember(bars_all_slivers,bars,'rows');
        ind_sliver_bars = ind(ind > 0);
        [~,ind_temp]    = unique(ind_sliver_bars);
        repeat_ind      = setdiff(1:length(ind_sliver_bars),ind_temp);
        repeat_values   = ind_sliver_bars(repeat_ind,1);
        shared_bars     = bars(unique(repeat_values,'stable'),:); % shared bars for this sliver
        num_shared_bars = size(shared_bars,1);   % number of shared bars for this sliver        
        
        %==============================================================================================================
        % COMPUTE POSSIBLE TRIANGLES AND THEIR QUALITY FACTORS
        %==============================================================================================================
        % CASE 1: Triangles associated to 1 vertex and 2 midpoints, e.g. triangle formed by nodes 1, midpoint in 2-3 and modipoint in 3-4
        %
        %               4                                                   4                      
        %               .                                                   .                      
        %            .  |     .                                          .  |     .      midpoint 3-4           
        %         .     |          .                                  .    _|_____----x           
        %      .         |             .                           . _ ----  |        |   .       
        %   1. - - - - - - - - - - - - - -. 3  ------------>    1x - - - - - - - - - - - - - -. 3  
        %      .          |            .                           . -___     |       |    .       
        %         .       |         .                                 .  ----_|___   |  .          
        %            .     |     .                                       .     |  ---x             
        %               .  |  .                                             .  |  .    midpoint 2-3            
        %                  .                                                   .                   
        %                  2                                                   2                   
        %
        [GCOORD_tri_1,q1] = tri_1_vertex_2_midpoints(GCOORD,barmid,sorted_nodes,I_fix);
        
        % CASE 2: Triangles associated to 3 vertex (tetrahedron faces), e.g., triangle formed by nodes 1, 2 and 3
        %
        %               4                                                   4                      
        %               .                                                   .                      
        %            .  |     .                                          .  |     .                 
        %         .     |          .                                  .     |          .           
        %      .         |             .                           .         |             .       
        %   1. - - - - - - - - - - - - - -. 3  ------------>    1x - - - - - - - - - - - - - -x 3  
        %      .          |            .                           .          |            .       
        %         .       |         .                                 .       |         .          
        %            .     |     .                                       .     |     .             
        %               .  |  .                                             .  |  .                         
        %                  .                                                   x                   
        %                  2                                                   2                   
        %
        q2 = tri_3_vertices(GCOORD,bars,I_fix);
        
        % CASE 3: Triangle associated to 2 vertices and 1 midpoint (in the shared bar), e.g., if the shared edge is 1-2, 
        % the triangle is formed by midpoint in 1-2, node 3 and node 4
        %
        %               4                                                   4                      
        %               .                                                   x                      
        %            .  |     .                                          . ||     .                 
        %         .     |          .                                  .   | |          .           
        %      .         |             .                           .      |  |             .       
        %   1. - - - - - - - - - - - - - -. 3  ------------>    1. - - - - - - - - - - - - - -x 3  
        %      .          |            .                           .     |   |     ___--- .       
        %         .       |         .                                 .  |   _|_---     .          
        %            .     |     .                                       x---  |     .             
        %               .  |  .                             midpoint 1-2    .  |  .                         
        %                  .                                                   .                   
        %                  2                                                   2                   
        %
        if num_shared_bars == 1 && (num_bnd_nodes == 0 || num_bnd_nodes == 1)% 1 shared edge and either no bnd nodes or 1 bnd node
            % Triangle associated to 2 vertices and 1 midpoint (of the shared bar)
            q3 = tri_2_vertices_1_midpoint(GCOORD,barmid,sorted_nodes,bars,shared_bars,I_fix);
        else
            q3 = 0;
        end
        
        %====================================================================================================================================
        % COMPARE TRIANGLES, CHOOSE THE BEST ONE, COMPUTE NODES TO BE ADDED AND REMOVED DEPENDING ON THE TYPE OF NODE (BOUNDARY OR INTERIOR)
        %====================================================================================================================================
        if num_bnd_nodes == 0 % interior sliver
            % Compare the quality factor of these 8 or 9 triangles
            [int_nodes_to_remove,pint_new] = ...
                compare_q_slivers_int_nodes(GCOORD_tri_1,q1,q2,q3,sorted_nodes,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix);
        
        elseif num_bnd_nodes == 1 % boundary sliver with 1 bnd node
            % Compare the quality factor of these 8 or 9 triangles
            if ~isempty(I_bnd1) % sliver on boundary 1
                [bnd1_nodes_to_remove,int_nodes_to_remove,pint_new] = ...
                    compare_q_slivers_1bnd_node(GCOORD_tri_1,q1,q2,q3,sorted_nodes,bnd1_nodes_to_remove,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix,I_bnd1,I_int);
            else                % sliver on boundary 2
                [bnd2_nodes_to_remove,int_nodes_to_remove,pint_new] = ...
                    compare_q_slivers_1bnd_node(GCOORD_tri_1,q1,q2,q3,sorted_nodes,bnd2_nodes_to_remove,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix,I_bnd2,I_int);
            end
            
        elseif num_bnd_nodes == 2 && num_shared_bars == 0 % boundary sliver with 2 bnd nodes but no shared bars
            % Compare the quality factor of these 8 triangles
            if ~isempty(I_bnd1) % sliver on boundary 1
                [bnd1_nodes_to_remove,int_nodes_to_remove,pbnd1_new,pint_new] = ...
                    compare_q_slivers_2bnd_nodes(GCOORD,GCOORD_tri_1,q1,q2,sorted_nodes,bnd1_nodes_to_remove,int_nodes_to_remove,pbnd1_new,pint_new,I_fix,I_bnd1,I_int);
            else                % sliver on boundary 2
                [bnd2_nodes_to_remove,int_nodes_to_remove,pbnd2_new,pint_new] = ....
                    compare_q_slivers_2bnd_nodes(GCOORD,GCOORD_tri_1,q1,q2,sorted_nodes,bnd2_nodes_to_remove,int_nodes_to_remove,pbnd2_new,pint_new,I_fix,I_bnd2,I_int);
            end
            
        elseif num_bnd_nodes == 3 && num_shared_bars == 0 % boundary sliver with 3 bnd nodes (special case)
            % Being the sliver with nodes 1,2,3 and 4:
            % If e.g. 1, 2 and 4 are the boundary nodes, there is 1 possible triangle associated to 3 nodes (bnd nodes):
            % Triangle formed by nodes 1, 2 and 4
            %                                                                          bnd node  
            %               4                                                        4
            %               .                                                        x
            %            .  |     .                                               .  |     .
            %         .     |          .                                       .     |          .
            %      .         |             .                    bnd node    .         |             .
            %   1. - - - - - - - - - - - - - -. 3  ------------>         1x - - - - - - - - - - - - - -. 3
            %      .          |            .                                .          |            .
            %         .       |         .                                      .       |         .
            %            .     |     .                                            .     |     .
            %               .  |  .                                                  .  |  .
            %                  .                                                        x
            %                  2                                                        2
            %                                                                             bnd node
            %
            nodes_to_remove     = sorted_nodes(I_int)'; % remove the int node
            int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove(:)];
            int_nodes_to_remove = unique(int_nodes_to_remove);
        end
    else % (some node has been removed)
        % Do nothing more and go to the following sliver in the loop
    end
end

%==========================================================================
% ADDING OR REJECTING NODES
%==========================================================================
% BOUNDARY 1
pbnd1(bnd1_nodes_to_remove-nfix,:) = []; % Remove bnd1 nodes
if ~isempty(pbnd1_new) && SETTINGS.r_int > 0
    % Project these new bnd1 nodes onto the spherical boundary1
    pbnd1_new_SPH      = cartesian2spherical(pbnd1_new);               % convert to spherical coordinates
    pbnd1_new_SPH(:,3) = SETTINGS.r_int*ones(size(pbnd1_new_SPH,1),1); % change the r-coordinate by r_int (projecting onto bnd 1)
    pbnd1_new          = spherical2cartesian(pbnd1_new_SPH);           % Cartesian coordinates of new bnd1 nodes after projecting on bnd1
end
pbnd1 = [pbnd1; pbnd1_new]; % Add new bnd1 nodes

% BOUNDARY 2
pbnd2(bnd2_nodes_to_remove-nfix-nbnd1,:) = []; % Remove bnd2 nodes
if ~isempty(pbnd2_new) && ~isempty(SETTINGS.r_ext)
    % Project these new bnd2 nodes onto the spherical boundary2
    pbnd2_new_SPH      = cartesian2spherical(pbnd2_new);               % convert to spherical coordinates
    pbnd2_new_SPH(:,3) = SETTINGS.r_ext*ones(size(pbnd2_new_SPH,1),1); % change the r-coordinate by r_ext (projecting onto bnd 2)
    pbnd2_new          = spherical2cartesian(pbnd2_new_SPH);           % Cartesian coordinates of new bnd2 nodes after projecting on bnd2
end
pbnd2 = [pbnd2; pbnd2_new]; % Add new bnd2 nodes

% INTERIOR
pint(int_nodes_to_remove-nfix-nbnd1-nbnd2,:) = [];               % Remove interior nodes
pint                                         = [pint; pint_new]; % Add new interior points

%==========================================================================
% DATA FOR PLOTS
%==========================================================================
GCOORD_new       = [pfix; pbnd1; pbnd2; pint];
EL2NOD_new       = delaunay(GCOORD_new);
% Remove elements created inside the interior boundary (boundary 1)
if SETTINGS.r_int > 0
    GCOORD_SPH_new   = cartesian2spherical(GCOORD_new);
    EL2NOD_new       = EL2NOD_new(~(sum(ismember(EL2NOD_new,find(abs(GCOORD_SPH_new(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
end
q_new            = tetra_mesh_quality(GCOORD_new,EL2NOD_new);
new_slivers      = EL2NOD_new(q_new < q_sliver,:);
sliver_nodes     = SLIVER2NOD_sorted(:);
sliver_bnd_nodes = sliver_nodes(sliver_nodes<=nfix+nbnd1+nbnd2);

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(44)
    hold on
    if nargin ~= 0
        clf
        plot_sphere(SETTINGS.r_ext)
    end
    faceColor = [0.6875 0.8750 0.8984];
    tetramesh(SLIVER2NOD_sorted,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
    if ~isempty(pbnd1_new)
        scatter3(pbnd1_new(:,1),pbnd1_new(:,2),pbnd1_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) % bnd1 node added in blue
    end
    if ~isempty(pbnd2_new)
        scatter3(pbnd2_new(:,1),pbnd2_new(:,2),pbnd2_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd2 node added in yellow
    end
    if ~isempty(pint_new)
        scatter3(pint_new(:,1),pint_new(:,2),pint_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]) % int node added in green
    end
    scatter3(GCOORD(sliver_bnd_nodes,1),GCOORD(sliver_bnd_nodes,2),GCOORD(sliver_bnd_nodes,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1])     % sliver bnd nodes in cyan
    if SETTINGS.r_int > 0
        scatter3(GCOORD(bnd1_nodes_to_remove,1),GCOORD(bnd1_nodes_to_remove,2),GCOORD(bnd1_nodes_to_remove,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]) % bnd1 node removed in white
    end
    scatter3(GCOORD(bnd2_nodes_to_remove,1),GCOORD(bnd2_nodes_to_remove,2),GCOORD(bnd2_nodes_to_remove,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1])     % bnd2 node removed in magenta
    scatter3(GCOORD(int_nodes_to_remove,1),GCOORD(int_nodes_to_remove,2),GCOORD(int_nodes_to_remove,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])     % int node removeds in red
    if SETTINGS.r_int > 0
        title({'MESH BEFORE FIX SLIVER (slivers in blue);'; ...
            ' bnd1 nodes: added (blue), removed(white); bnd2 nodes: added (yellow), removed (magenta); int nodes: added (green), removed (red)'})
    else
        title({'MESH BEFORE FIX SLIVER (slivers in blue);'; ...
            ' bnd2 nodes: added (yellow), removed (magenta); int nodes: added (green), removed (red)'})
    end
    if nargin ~= 0
        text(6000,-6000,6000,['worst q = ',sprintf('%.3f',min(q))], ...
            'HorizontalAlignment','center','BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black')
    else
        text(1,0,1.2,['worst q = ',sprintf('%.3f',min(q))], ...
            'HorizontalAlignment','center','BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black')
    end
    
    figure(45)
    clf
    hold on
    if nargin ~= 0
        plot_sphere(SETTINGS.r_ext)
    else
        axis equal
        view(142.5,30)
        xlabel('X');ylabel('Y');zlabel('Z')
        scatter3(GCOORD_new(EL2NOD_new,1),GCOORD_new(EL2NOD_new,2),GCOORD_new(EL2NOD_new,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
        tetramesh(EL2NOD_new,GCOORD_new,'FaceColor',faceColor,'FaceAlpha',0);
    end
    tetramesh(new_slivers,GCOORD_new,'FaceColor',faceColor,'FaceAlpha',0.3);
    if ~isempty(pbnd1_new)
        scatter3(pbnd1_new(:,1),pbnd1_new(:,2),pbnd1_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) % bnd1 node added in blue
    end
    if ~isempty(pbnd2_new)
        scatter3(pbnd2_new(:,1),pbnd2_new(:,2),pbnd2_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd2 node added in yellow
    end
    if ~isempty(pint_new)
        scatter3(pint_new(:,1),pint_new(:,2),pint_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]) % int node added in green
    end
    scatter3(GCOORD(sliver_bnd_nodes,1),GCOORD(sliver_bnd_nodes,2),GCOORD(sliver_bnd_nodes,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1])     % sliver bnd nodes in cyan
    if SETTINGS.r_int > 0
        scatter3(GCOORD(bnd1_nodes_to_remove,1),GCOORD(bnd1_nodes_to_remove,2),GCOORD(bnd1_nodes_to_remove,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]) % bnd1 node removed in white
    end
    scatter3(GCOORD(bnd2_nodes_to_remove,1),GCOORD(bnd2_nodes_to_remove,2),GCOORD(bnd2_nodes_to_remove,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1])     % bnd2 node removed in magenta
    scatter3(GCOORD(int_nodes_to_remove,1),GCOORD(int_nodes_to_remove,2),GCOORD(int_nodes_to_remove,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])     % int node removed in red
    if SETTINGS.r_int > 0
        title({'MESH AFTER FIX SLIVER (slivers in blue);'; ...
            ' bnd1 nodes: added (blue), removed(white); bnd2 nodes: added (yellow), removed (magenta); int nodes: added (green), removed (red)'})
    else
        title({'MESH AFTER FIX SLIVER (slivers in blue);'; ...
            ' bnd2 nodes: added (yellow), removed (magenta); int nodes: added (green), removed (red)'})
    end
    if nargin ~= 0
        text(6000,-6000,6000,['worst q = ',sprintf('%.3f',min(q_new))], ...
            'HorizontalAlignment','center','BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black')
    else
        text(1,0,1.2,['worst q = ',sprintf('%.3f',min(q_new))], ...
            'HorizontalAlignment','center','BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black')
    end
end

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd1                     = pbnd1;
MESH.pbnd2                     = pbnd2;
MESH.pint                      = pint;
MESH.GCOORD                    = GCOORD_new;
MESH.EL2NOD                    = EL2NOD_new;
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
MESH.rel_change                = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
MESH.rel_change_abs            = abs(MESH.rel_change);        % absolute value of relative bar-length change
MESH.mean_misfit_bar_length    = sum(MESH.rel_change_abs)/size(MESH.L,1);
MESH.q                         = q_new;
MESH.s                         = shape_measure(MESH.GCOORD,MESH.EL2NOD);

end % END OF FUNCTION fix_slivers

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [GCOORD_tri,q_tri] = tri_1_vertex_2_midpoints(GCOORD,barmid,sorted_nodes,I_fix)
% Usage: [GCOORD_tri,q_tri] = tri_1_vertex_2_midpoints(GCOORD,barmid,sorted_nodes,I_fix)

% Purpose:
%   Given a sliver, compute 4 triangles associated to 1 node and 2
%   midpoints, and their quality factors. For example, being the sliver
%   with nodes 1,2,3 and 4, this routine compute the following triangles:
%   
%               4
%               .                        
%            .  |     .                  
%         .     |          .            
%      .         |             .        
%   1. - - - - - - - - - - - - - -. 3   
%      .          |            .
%         .       |         .       
%            .     |     .             
%               .  |  .
%                  .                  
%                  2
%   
%   Triangle formed by nodes 1, midpoint in 2-3 and modipoint in 3-4
%   Triangle formed by nodes 2, midpoint in 1-4 and modipoint in 3-4
%   Triangle formed by nodes 3, midpoint in 1-2 and modipoint in 1-4
%   Triangle formed by nodes 4, midpoint in 1-2 and modipoint in 2-3
%
%   NOTE: Fixed nodes cannot be removed.
%   If there is 1 fixed node, e.g. node 2, there is only one possible 
%   triangle, formed by nodes 2, midpoint in 1-4 and modipoint in 3-4
%
% Input:
%   GCOORD       : [matrix] : Cartesian coordinates of the mesh
%   barmid       : [matrix] : coordiantes of the midpoints of each bar in
%                             the sliver
%   sorted_nodes : [vector] : sliver nodes (vertices) sorted from the
%                             smallest one
%   I_fix        : [vector] : local index of bnd nodes for this sliver
%
% Output:
%   GCOORD_tri   : [matrix] : coordinates of the 4 triangles
%   q_tri        : [vector] : quality factor for each triangle
%
% JMT Jul 2016

GCOORD_tri  = zeros(12,3);
q_tri  = zeros(4,1);
for j = 1:4
    % distance from the node j in the sliver to the 6 midpoints of the bars (each tetraedron has 6 bars)
    d        = sqrt(sum((barmid-repmat(GCOORD(sorted_nodes(j),:),6,1)).^2,2));
    [~,id]   = sort(d,'descend'); % sort the distance in order to get the 2 largest ones
    GCOORD_tri(3*j-2:3*j,:) = [GCOORD(sorted_nodes(j),:); barmid(id(1:2),:)]; % triangle associated to node j
    a        = d(id(1));
    b        = d(id(2));
    c        = sqrt(sum((barmid(id(2),:)-barmid(id(1),:)).^2,2));
    q_tri(j) = ((b + c - a)*(c + a - b)*(a + b - c)) / (a*b*c); % quality factor for the triangle associated to node j
end

if size(I_fix,2) == 1 % if there is 1 fixed node, select the triangle formed by the fixed node and the opposite midpoints
    GCOORD_tri = GCOORD_tri(3*I_fix-2:3*I_fix,:); % triangle associated to node I_fix
    q_tri      = q_tri(I_fix); % quality factor for the triangle associated to node I_fix
elseif size(I_fix,2) > 1 % if there is 2 or more fixed nodes, do nothing, since this routine requires removing 3 nodes
    GCOORD_tri = [];
    q_tri      = 0;
end

end % END OF SUBFUNCTION tri_1_node_2_midpoints

function [q_tri] = tri_3_vertices(GCOORD,bars,I_fix)
% Usage: [q_tri] = tri_3_vertices(GCOORD,bars,I_fix)

% Purpose:
%   Given a sliver, compute 4 triangles associated to 3 nodes of each face,
%   and their quality factors. For example, being the sliver with nodes
%   1,2,3 and 4, this routine compute the following triangles:
%   
%               4
%               .                        
%            .  |     .                  
%         .     |          .            
%      .         |             .        
%   1. - - - - - - - - - - - - - -. 3   
%      .          |            .
%         .       |         .       
%            .     |     .             
%               .  |  .
%                  .                  
%                  2
%   
%   Triangle formed by nodes 1, 2 and 3
%   Triangle formed by nodes 1, 2 and 4
%   Triangle formed by nodes 1, 3 and 4
%   Triangle formed by nodes 2, 3 and 4
%
% Input:
%   GCOORD : [matrix] : Cartesian coordinates of the mesh
%   bars   : [matrix] : bars list
%   I_fix  : [vector] : local index of bnd nodes for this sliver
%
% Output:
%   q_tri  : [vector] : quality factor for each triangle
%
% JMT Jul 2016

d_bars = sqrt(sum((GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:)).^2,2));
%         face 123   face 124   face 134   face 234
a      = [d_bars(1); d_bars(1); d_bars(2); d_bars(4)];
b      = [d_bars(4); d_bars(5); d_bars(6); d_bars(6)];
c      = [d_bars(2); d_bars(3); d_bars(3); d_bars(5)];
q_tri  = ((b + c - a).*(c + a - b).*(a + b - c)) ./ (a.*b.*c);

if size(I_fix,2) == 1 % if there is 1 fixed node, select the triangles that contain the fixed node
    if I_fix == 1 % local node 1 is fixed -> triangles 123, 124 and 134
        q_tri(4) = []; % remove triangle 234 (since it would imply removing node 1 which is fixed and cannot be removed)
    elseif I_fix == 2 % local node 2 is fixed -> triangles 123, 124 and 234
        q_tri(3) = []; % remove triangle 134 (since it would imply removing node 2 which is fixed and cannot be removed)
    elseif I_fix == 3 % local node 3 is fixed -> triangles 123, 134 and 234
        q_tri(2) = []; % remove triangle 124 (since it would imply removing node 3 which is fixed and cannot be removed)
    elseif I_fix == 4 % local node 4 is fixed -> triangles 124, 134 and 234
        q_tri(1) = []; % remove triangle 123 (since it would imply removing node 4 which is fixed and cannot be removed)
    end
elseif size(I_fix,2) == 2 % if there is 2 fixed nodes, select the triangles that contain the fixed nodes
    if sum(ismember([1 2],I_fix)) == 2
        q_tri = q_tri([1 2]);
    elseif sum(ismember([1 3],I_fix)) == 2
        q_tri = q_tri([1 3]);
    elseif sum(ismember([1 4],I_fix)) == 2
        q_tri = q_tri([2 3]);
    elseif sum(ismember([2 3],I_fix)) == 2
        q_tri = q_tri([1 4]);
    elseif sum(ismember([2 4],I_fix)) == 2
        q_tri = q_tri([2 4]);
    elseif sum(ismember([3 4],I_fix)) == 2
        q_tri = q_tri([3 4]);
    end
elseif size(I_fix,2) == 3 % if there is 3 fixed nodes, select the triangle that contains the fixed nodes
    if sum(ismember([1 2 3],I_fix)) == 3
        q_tri = q_tri(1);
    elseif sum(ismember([1 2 4],I_fix)) == 3
        q_tri = q_tri(2);
    elseif sum(ismember([1 3 4],I_fix)) == 3
        q_tri = q_tri(3);
    elseif sum(ismember([2 3 4],I_fix)) == 3
        q_tri = q_tri(4);
    end
elseif size(I_fix,2) == 4 % if there is 4 fixed nodes, do nothing, since this routine requires removing 1 node
    q_tri      = 0;
end

end % END OF SUBFUNCTION tri_3_nodes

function [q_tri] = tri_2_vertices_1_midpoint(GCOORD,barmid,sorted_nodes,bars,shared_bars,I_fix)
% Usage: [q_tri] = tri_2_vertices_1_midpoint(GCOORD,barmid,bars,sorted_nodes,shared_bars,I_fix)

% Purpose:
%   Given a sliver which is sharing a bar (spring) with anoter sliver,
%   compute the triangle associated to the midpoint in the shared bar and  
%   the opposite 2 nodes, and its quality factor. For example, being the 
%   sliver with nodes 1,2,3 and 4  and sharing the edge 1-2, this routine 
%   compute the following triangle: 
%   
%               4
%               .                        
%            .  |     .                  
%         .     |          .            
%      .         |             .        
%   1. - - - - - - - - - - - - - -. 3   
%      .          |            .
%         .       |         .       
%            .     |     .             
%               .  |  .
%                  .                  
%                  2
%   
%   Triangle formed by midpoint between 1-2, node 3 and node 4 
%
% Input:
%   GCOORD       : [matrix] : Cartesian coordinates of the mesh
%   barmid       : [matrix] : coordiantes of the midpoints of each bar in
%                             the sliver
%   sorted_nodes : [vector] : sliver nodes (vertices) sorted from the
%                             smallest one
%   bars         : [matrix] : bars list
%   shared_bars  : [matrix] : shared bars list
%   I_fix        : [vector] : local index of bnd nodes for this sliver
%
% Output:
%   q_tri        : [vector] : quality factor for each triangle
%
% JMT Jul 2016

[non_shared_nodes,~] = setdiff(sorted_nodes,shared_bars);
a      = sqrt(sum((GCOORD(non_shared_nodes(1),:)-barmid(ismember(bars,shared_bars,'rows'),:)).^2,2)); % distance from midpoint of shared bar to one of the non shared nodes
b      = sqrt(sum((GCOORD(non_shared_nodes(2),:)-barmid(ismember(bars,shared_bars,'rows'),:)).^2,2)); % distance from midpoint of shared bar to the other non shared node
c      = sqrt(sum((GCOORD(non_shared_nodes(2),:)-GCOORD(non_shared_nodes(1),:)).^2,2));               % distance between the non shared nodes
q_tri  = ((b + c - a)*(c + a - b)*(a + b - c)) / (a*b*c);

if size(I_fix,2) > 0 % if there is 1 or more fixed nodes, do nothing, since this routine requires removing 2 nodes that share a bar
    q_tri      = 0;
end

end % END OF SUBFUNCTION tri_2_vertices_1_midpoint

function [int_nodes_to_remove,pint_new] = ...
    compare_q_slivers_int_nodes(GCOORD_tri_1,q1,q2,q3,sorted_nodes,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix)
% Usage: [int_nodes_to_remove,pint_new] = ...
%   compare_q_slivers_int_nodes(GCOORD_tri_1,q1,q2,q3,sorted_nodes,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix)
%
% Purpose:
%   Compare quality factor of possible triangles for a sliver with interior
%   nodes. Then choose the best triangle and select which nodes have to be
%   removed and create nodes to be added.
%
% Input:
%   GCOORD_tri_1        : [matrix] : coordinates of the 4 triangles for the
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q1                  : [vector] : quality factor of each triangle for
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q2                  : [vector] : quality factor of each triangle for
%                                    case 2 (triangles associated to 3
%                                    nodes of each face)
%   q3                  : [scalar] : quality factor of the triangle for 
%                                    case 3 (triangle associated to the
%                                    midpoint in the shared bar and the
%                                    opposite 2 nodes)
%   sorted_nodes        : [vector] : sorted nodes for this sliver from
%                                    smallest one to the biggest one 
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pint_new            : [matrix] : coordinates for new interior nodes
%   shared_bars         : [matrix] : shared bars for this sliver
%   barmid              : [matrix] : midpoints of the bars for this sliver
%   bars                : [matrix] : bars for this sliver
%   I_fix               : [vector] : local index of bnd nodes for this
%                                    sliver
%
% Output:
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pint_new            : [matrix] : coordinates for new interior nodes
%
% JMT Jul 2016

if max(q1) > max([q3; q2]) % -> triangle associated to 1 node and 2 midpoints
    if isempty(I_fix)
        [~,iq_tri]                 = sort(q1,'descend'); % sort quality factor in order to get the best triangle and the associated node --> iq_tri(1)
        nodes_to_remove            = sorted_nodes';
        nodes_to_remove(iq_tri(1)) = []; % remove from this list the node that we want to keep --> iq_tri(1)
        int_nodes_to_remove        = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove        = unique(int_nodes_to_remove(:));
        p_to_add                   = GCOORD_tri_1(3*iq_tri(1)-1:3*iq_tri(1),:); % midpoints for the 2 bars opposite to the node iq_tri(1)
        pint_new                   = [pint_new; p_to_add];
        pint_new                   = unique(pint_new,'rows','stable');
    else
        nodes_to_remove        = sorted_nodes';
        nodes_to_remove(I_fix) = []; % remove from this list the node that we want to keep --> I_fix
        int_nodes_to_remove    = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove    = unique(int_nodes_to_remove(:));
        p_to_add               = GCOORD_tri_1(2:3,:); % midpoints for the 2 bars opposite to the node I_fix
        pint_new               = [pint_new; p_to_add];
        pint_new               = unique(pint_new,'rows','stable');
    end
    
elseif max(q2) > max([q1; q3]) % -> triangle associated to 3 nodes of each face
    [~,iq_tri2] = sort(q2,'descend'); % sort the quality factor in order to get the best triangle and the associated node
    if isempty(I_fix)
        if     iq_tri2(1) == 1 % -> face 123
            nodes_to_remove = sorted_nodes(4); % remove local node 4
        elseif iq_tri2(1) == 2 % -> face 124
            nodes_to_remove = sorted_nodes(3); % remove local node 3
        elseif iq_tri2(1) == 3 % -> face 134
            nodes_to_remove = sorted_nodes(2); % remove local node 2
        elseif iq_tri2(1) == 4 % -> face 234
            nodes_to_remove = sorted_nodes(1); % remove local node 1
        end
        int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove = unique(int_nodes_to_remove(:));
    elseif size(I_fix,2) == 1 % there is 1 fixed node
        if I_fix == 1 % local node 1 is fixed -> triangles 123, 124 and 134
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 3 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            end
        elseif I_fix == 2 % local node 2 is fixed -> triangles 123, 124 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif I_fix == 3 % local node 3 is fixed -> triangles 123, 134 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif I_fix == 4 % local node 4 is fixed -> triangles 124, 134 and 234
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        end
        int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove = unique(int_nodes_to_remove(:));
    elseif size(I_fix,2) == 2 % there are 2 fixed nodes
        if sum(ismember([1 2],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            end
        elseif sum(ismember([1 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            end
        elseif sum(ismember([1 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            end
        elseif sum(ismember([2 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif sum(ismember([2 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif sum(ismember([3 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        end
        int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove = unique(int_nodes_to_remove(:));
    elseif size(I_fix,2) == 3 % there are 3 fixed nodes
        if sum(ismember([1 2 3],I_fix)) == 3
            nodes_to_remove = sorted_nodes(4); % remove local node 4
        elseif sum(ismember([1 2 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(3); % remove local node 3
        elseif sum(ismember([1 3 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(2); % remove local node 2
        elseif sum(ismember([2 3 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(1); % remove local node 1
        end
        int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove = unique(int_nodes_to_remove(:));
    elseif size(I_fix,2) == 4 % there are 4 fixed nodes (do nothing)
    end
    
elseif max(q3) > max([q1; q2]) % -> triangle associated to 1 midpoint in the shared bar and the opposite 2 nodes
    nodes_to_remove     = shared_bars'; % remove nodes in the shared bar
    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
    int_nodes_to_remove = unique(int_nodes_to_remove(:));
    p_to_add            = barmid(ismember(bars,shared_bars,'rows'),:); % add midpoint in the shared bar
    pint_new            = [pint_new; p_to_add];
    pint_new            = unique(pint_new,'rows','stable');
end

end % END OF SUBFUNCTION compare_q_slivers_int_nodes

function [bnd_nodes_to_remove,int_nodes_to_remove,pint_new] = ...
    compare_q_slivers_1bnd_node(GCOORD_tri_1,q1,q2,q3,sorted_nodes,bnd_nodes_to_remove,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix,I_bnd,I_int)
% Usage: [bnd_nodes_to_remove,int_nodes_to_remove,pint_new] = ...
%   compare_q_slivers_1bnd_node(GCOORD_tri_1,q1,q2,q3,sorted_nodes,bnd_nodes_to_remove,int_nodes_to_remove,pint_new,shared_bars,barmid,bars,I_fix,I_bnd,I_int)
%
% Purpose:
%   Compare quality factor of possible triangles for a sliver with 1
%   boundary node. Then choose the best triangle and select which nodes
%   have to be removed and create nodes to be added.
%
% Input:
%   GCOORD_tri_1        : [matrix] : coordinates of the 4 triangles for the
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q1                  : [vector] : quality factor of each triangle for
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q2                  : [vector] : quality factor of each triangle for
%                                    case 2 (triangles associated to 3
%                                    nodes of each face)
%   q3                  : [scalar] : quality factor of the triangle for 
%                                    case 3 (triangle associated to the
%                                    midpoint in the shared bar and the
%                                    opposite 2 nodes)
%   sorted_nodes        : [vector] : sorted nodes for this sliver from
%                                    smallest one to the biggest one
%   bnd_nodes_to_remove : [vector] : list of bnd nodes to remove
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pint_new            : [matrix] : coordinates for new interior nodes
%   shared_bars         : [matrix] : shared bars for this sliver
%   barmid              : [matrix] : midpoints of the bars for this sliver
%   bars                : [matrix] : bars for this sliver
%   I_fix               : [vector] : local index of fix nodes for this
%                                    sliver
%   I_bnd               : [vector] : local index of bnd nodes for this
%                                    sliver
%   I_int               : [vector] : local index of int nodes for this
%                                    sliver
%
% Output:
%   bnd_nodes_to_remove : [vector] : list of bnd nodes to remove
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pint_new            : [matrix] : coordinates for new interior nodes
%
% JMT Jul 2016

if max(q1) > max([q3; q2]) % -> triangle associated to 1 node and 2 midpoints
    [~,iq_tri] = sort(q1,'descend'); % sort quality factor in order to get the best triangle and the associated node --> iq_tri(1)
    if ismember(iq_tri(1),I_int) % means that the node we want to keep (node associated to the best triangle) is an interior node
        % remove bnd node
        temp_bnd_nodes_to_remove = sorted_nodes(I_bnd)'; % bnd node in the sliver
        bnd_nodes_to_remove      = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove];
        bnd_nodes_to_remove      = unique(bnd_nodes_to_remove);
        % remove interior nodes
        temp_nodes_to_remove            = sorted_nodes';
        temp_nodes_to_remove(iq_tri(1)) = []; % remove from this list the node we want to keep (node associated to the best triangle)
        temp_int_nodes_to_remove        = temp_nodes_to_remove(~ismember(temp_nodes_to_remove,temp_bnd_nodes_to_remove)); % remove the bnd node from this list
        int_nodes_to_remove             = [int_nodes_to_remove; temp_int_nodes_to_remove];
        int_nodes_to_remove             = unique(int_nodes_to_remove(:));
        % add interior nodes
        p_to_add                        = GCOORD_tri_1(3*iq_tri(1)-1:3*iq_tri(1),:); % midpoints for the 2 bars opposite to the node iq_tri(1)
        pint_new                        = [pint_new; p_to_add];
        pint_new                        = unique(pint_new,'rows','stable');
    elseif ismember(iq_tri(1),I_bnd) % means that the node we want to keep (node associated to the best triangle) is a bnd node
        % remove interior nodes
        temp_int_nodes_to_remove            = sorted_nodes';
        temp_int_nodes_to_remove(iq_tri(1)) = []; % remove from this list the node that we want to keep --> iq_tri(1)
        int_nodes_to_remove                 = [int_nodes_to_remove; temp_int_nodes_to_remove];
        int_nodes_to_remove                 = unique(int_nodes_to_remove(:));
        % add interior nodes
        p_to_add                            = GCOORD_tri_1(3*iq_tri(1)-1:3*iq_tri(1),:); % midpoints for the 2 bars opposite to the node iq_tri(1)
        pint_new                            = [pint_new; p_to_add];
        pint_new                            = unique(pint_new,'rows','stable');
    elseif ismember(iq_tri(1),I_fix) % means that the node we want to keep (node associated to the best triangle) is a fixed node
        % remove bnd node
        temp_bnd_nodes_to_remove = sorted_nodes(I_bnd)'; % bnd node in the sliver
        bnd_nodes_to_remove      = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove];
        bnd_nodes_to_remove      = unique(bnd_nodes_to_remove);
        % remove interior nodes
        temp_nodes_to_remove            = sorted_nodes';
        temp_nodes_to_remove(iq_tri(1)) = []; % remove from this list the node we want to keep (node associated to the best triangle)
        temp_int_nodes_to_remove        = temp_nodes_to_remove(~ismember(temp_nodes_to_remove,temp_bnd_nodes_to_remove)); % remove the bnd node from this list
        int_nodes_to_remove             = [int_nodes_to_remove; temp_int_nodes_to_remove];
        int_nodes_to_remove             = unique(int_nodes_to_remove(:));
        % add interior nodes
        p_to_add                        = GCOORD_tri_1(2:3,:); % midpoints for the 2 bars opposite to the node I_fix
        pint_new                        = [pint_new; p_to_add];
        pint_new                        = unique(pint_new,'rows','stable');
    end
    
elseif max(q2) > max([q1; q3]) % -> triangle associated to 3 nodes of each face
    [~,iq_tri2] = sort(q2,'descend'); % sort the quality factor in order to get the best triangle and the associated node
    if isempty(I_fix)
        if     iq_tri2(1) == 1 % -> face 123
            nodes_to_remove = sorted_nodes(4); % remove local node 4 (check if local node 4 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif iq_tri2(1) == 2 % -> face 124
            nodes_to_remove = sorted_nodes(3); % remove local node 3 (check if local node 3 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif iq_tri2(1) == 3 % -> face 134
            nodes_to_remove = sorted_nodes(2); % remove local node 2 (check if local node 2 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif iq_tri2(1) == 4 % -> face 234
            nodes_to_remove = sorted_nodes(1); % remove local node 1 (check if local node 1 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        end
    elseif size(I_fix,2) == 1 % there is 1 fixed node
        if I_fix == 1 % local node 1 is fixed -> triangles 123, 124 and 134
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 2 % local node 2 is fixed -> triangles 123, 124 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 3 % local node 3 is fixed -> triangles 123, 134 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 4 % local node 4 is fixed -> triangles 124, 134 and 234
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        end
    elseif size(I_fix,2) == 2 % there are 2 fixed nodes
        if sum(ismember([1 2],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif sum(ismember([1 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif sum(ismember([1 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif sum(ismember([2 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif sum(ismember([2 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif sum(ismember([3 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if sorted_nodes(I_bnd) == nodes_to_remove
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        end
    elseif size(I_fix,2) == 3 % there are 3 fixed nodes    
        if sum(ismember([1 2 3],I_fix)) == 3
            nodes_to_remove = sorted_nodes(4); % remove local node 4
            % (check if local node 4 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif sum(ismember([1 2 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(3); % remove local node 3
            % (check if local node 3 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif sum(ismember([1 3 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(2); % remove local node 2
            % (check if local node 2 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        elseif sum(ismember([2 3 4],I_fix)) == 3
            nodes_to_remove = sorted_nodes(1); % remove local node 1
            % (check if local node 1 is bnd or int node)
            if sorted_nodes(I_bnd) == nodes_to_remove
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove(:));
            end
        end
    end
    
elseif max(q3) > max([q1; q2]) % -> triangle associated to 1 midpoint in the shared bar and the opposite 2 nodes
    if ismember(sorted_nodes(I_bnd),shared_bars) % bnd node is on the shared edge
        temp_bnd_nodes_to_remove = sorted_nodes(I_bnd); % remove bnd node (which is on the shared edge)
        bnd_nodes_to_remove      = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove];
        bnd_nodes_to_remove      = unique(bnd_nodes_to_remove(:));
        
        temp_int_nodes_to_remove = shared_bars'; % remove int node (which is on the shared edge)
        temp_int_nodes_to_remove = temp_int_nodes_to_remove(~ismember(temp_int_nodes_to_remove,temp_bnd_nodes_to_remove));
        int_nodes_to_remove      = [int_nodes_to_remove; temp_int_nodes_to_remove];
        int_nodes_to_remove      = unique(int_nodes_to_remove(:));
        
        p_to_add            = barmid(ismember(bars,shared_bars,'rows'),:); % add midpoint (interior node) in the shared bar
        pint_new            = [pint_new; p_to_add];
        pint_new            = unique(pint_new,'rows','stable');
    else % bnd node is not on the shared edge
        nodes_to_remove     = shared_bars'; % remove int nodes (which are on the shared edge)
        int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
        int_nodes_to_remove = unique(int_nodes_to_remove(:));
        p_to_add            = barmid(ismember(bars,shared_bars,'rows'),:); % add midpoint (interior node) in the shared bar
        pint_new            = [pint_new; p_to_add];
        pint_new            = unique(pint_new,'rows','stable');
    end
end

end % END OF SUBFUNCTION compare_q_slivers_1bnd_node

function [bnd_nodes_to_remove,int_nodes_to_remove,pbnd_new,pint_new] = ...
    compare_q_slivers_2bnd_nodes(GCOORD,GCOORD_tri_1,q1,q2,sorted_nodes,bnd_nodes_to_remove,int_nodes_to_remove,pbnd_new,pint_new,I_fix,I_bnd,I_int)
% Usage: [bnd_nodes_to_remove,int_nodes_to_remove,pbnd_new,pint_new] = ...
%   compare_q_slivers_2bnd_nodes(GCOORD,GCOORD_tri_1,q1,q2,sorted_nodes,bnd_nodes_to_remove,int_nodes_to_remove,pbnd_new,pint_new,I_fix,I_bnd,I_int)
%
% Purpose:
%   Compare quality factor of possible triangles for a sliver with 2
%   boundary node. Then choose the best triangle and select which nodes
%   have to be removed and create nodes to be added.
%
% Input:
%   GCOORD              : [matrix] : coordinates of the mesh
%   GCOORD_tri_1        : [matrix] : coordinates of the 4 triangles for the
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q1                  : [vector] : quality factor of each triangle for
%                                    case 1 (triangles associated to 1 node
%                                    and 2 midpoints)
%   q2                  : [vector] : quality factor of each triangle for
%                                    case 2 (triangles associated to 3
%                                    nodes of each face)
%   sorted_nodes        : [vector] : sorted nodes for this sliver from
%                                    smallest one to the biggest one
%   bnd_nodes_to_remove : [vector] : list of bnd nodes to remove
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pbnd_new            : [matrix] : coordinates for new boundary nodes
%   pint_new            : [matrix] : coordinates for new interior nodes
%   I_fix               : [vector] : local index of fix nodes for this
%                                    sliver
%   I_bnd               : [vector] : local index of bnd nodes for this
%                                    sliver 
%   I_int               : [vector] : local index of int nodes for this
%                                    sliver 
%
% Output:
%   bnd_nodes_to_remove : [vector] : list of bnd nodes to remove
%   int_nodes_to_remove : [vector] : list of int nodes to remove
%   pbnd_new            : [matrix] : coordinates for new boundary nodes
%   pint_new            : [matrix] : coordinates for new interior nodes
%
% JMT Jul 2016

if max(q1) > max(q2) % -> triangle associated to 1 node and 2 midpoints
    [~,iq_tri] = sort(q1,'descend'); % sort quality factor in order to get the best triangle and the associated node --> iq_tri(1)
    if ismember(iq_tri(1),I_int) % means that the node we want to keep (node associated to the best triangle) is an interior node
        % remove bnd nodes
        temp_bnd_nodes_to_remove            = sorted_nodes(I_bnd)'; % bnd node in the sliver
        bnd_nodes_to_remove                 = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove];
        bnd_nodes_to_remove                 = unique(bnd_nodes_to_remove);
        % remove interior node
        temp_int_nodes_to_remove            = sorted_nodes;
        temp_int_nodes_to_remove(iq_tri(1)) = []; % remove from this list the node we want to keep (node associated to the best triangle)
        temp_int_nodes_to_remove            = temp_int_nodes_to_remove(~ismember(temp_int_nodes_to_remove,temp_bnd_nodes_to_remove)); % remove the bnd node from this list
        int_nodes_to_remove                 = [int_nodes_to_remove; temp_int_nodes_to_remove(:)];
        int_nodes_to_remove                 = unique(int_nodes_to_remove);
        % add nodes (1 bnd and 1 int)
        p_to_add                            = GCOORD_tri_1(3*iq_tri(1)-1:3*iq_tri(1),:); % midpoints for the 2 bars opposite to the node iq_tri(1)
        % check which p_to_add is interior and which one is on the boundary
        pbnd_to_check   = (GCOORD(sorted_nodes(I_bnd(1)),:) + GCOORD(sorted_nodes(I_bnd(2)),:))/2; % midpoint between 2 bnd nodes
        pbnd_new_to_add = p_to_add(ismember(p_to_add,pbnd_to_check))';
        pbnd_new        = [pbnd_new; pbnd_new_to_add];
        pbnd_new        = unique(pbnd_new,'rows','stable');
        pint_new_to_add = p_to_add;
        pint_new_to_add(ismember(p_to_add,pbnd_to_check)) = []; % remove pbnd_new
        pint_new        = [pint_new; pint_new_to_add];
        pint_new        = unique(pint_new,'rows','stable');
    elseif ismember(iq_tri(1),I_bnd) % means that the node we want to keep (node associated to the best triangle) is a bnd node
        % remove the other bnd node
        temp_bnd_nodes_to_remove = sorted_nodes(I_bnd)'; % bnd node in the sliver
        temp_bnd_nodes_to_remove = temp_bnd_nodes_to_remove(~ismember(sorted_nodes(I_bnd)',sorted_nodes(iq_tri(1))));
        bnd_nodes_to_remove      = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove(:)];
        bnd_nodes_to_remove      = unique(bnd_nodes_to_remove);
        % remove int nodes
        temp_int_nodes_to_remove = sorted_nodes';
        temp_int_nodes_to_remove = temp_int_nodes_to_remove(~ismember(temp_int_nodes_to_remove,sorted_nodes(I_bnd)));
        int_nodes_to_remove      = [int_nodes_to_remove; temp_int_nodes_to_remove(:)];
        int_nodes_to_remove      = unique(int_nodes_to_remove);
        % add int nodes
        p_to_add                 = GCOORD_tri_1(3*iq_tri(1)-1:3*iq_tri(1),:); % midpoints (interior) for the 2 bars opposite to the node iq_tri(1)
        pint_new                 = [pint_new; p_to_add];
        pint_new                 = unique(pint_new,'rows','stable');
    elseif ismember(iq_tri(1),I_fix) % means that the node we want to keep (node associated to the best triangle) is a fixed node
        % remove bnd nodes
        temp_bnd_nodes_to_remove            = sorted_nodes(I_bnd)'; % bnd node in the sliver
        bnd_nodes_to_remove                 = [bnd_nodes_to_remove; temp_bnd_nodes_to_remove];
        bnd_nodes_to_remove                 = unique(bnd_nodes_to_remove);
        % remove interior node
        temp_int_nodes_to_remove            = sorted_nodes;
        temp_int_nodes_to_remove(iq_tri(1)) = []; % remove from this list the node we want to keep (node associated to the best triangle)
        temp_int_nodes_to_remove            = temp_int_nodes_to_remove(~ismember(temp_int_nodes_to_remove,temp_bnd_nodes_to_remove)); % remove the bnd node from this list
        int_nodes_to_remove                 = [int_nodes_to_remove; temp_int_nodes_to_remove(:)];
        int_nodes_to_remove                 = unique(int_nodes_to_remove);
        % add nodes (1 bnd and 1 int)
        p_to_add                            = GCOORD_tri_1(2:3,:); % midpoints for the 2 bars opposite to the node I_fix
        % check which p_to_add is interior and which one is on the boundary
        pbnd_to_check   = (GCOORD(sorted_nodes(I_bnd(1)),:) + GCOORD(sorted_nodes(I_bnd(2)),:))/2; % midpoint between 2 bnd nodes
        pbnd_new_to_add = p_to_add(ismember(p_to_add,pbnd_to_check))';
        pbnd_new        = [pbnd_new; pbnd_new_to_add];
        pbnd_new        = unique(pbnd_new,'rows','stable');
        pint_new_to_add = p_to_add;
        pint_new_to_add(ismember(p_to_add,pbnd_to_check)) = []; % remove pbnd_new
        pint_new        = [pint_new; pint_new_to_add];
        pint_new        = unique(pint_new,'rows','stable');
    end
    
else % means that max(q_tri2) >= max(q_tri) -> triangle associated to 3 nodes of each face 
    [~,iq_tri2] = sort(q2,'descend'); % sort the quality factor in order to get the best triangle and the associated node
    if isempty(I_fix)
        if     iq_tri2(1) == 1 % -> face 123
            nodes_to_remove = sorted_nodes(4); % remove local node 4 (check if local node 4 is bnd or int node)
            if any(sorted_nodes(I_bnd) == nodes_to_remove)
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove);
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove);
            end
        elseif iq_tri2(1) == 2 % -> face 124
            nodes_to_remove = sorted_nodes(3); % remove local node 3 (check if local node 3 is bnd or int node)
            if any(sorted_nodes(I_bnd) == nodes_to_remove)
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove);
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove);
            end
        elseif iq_tri2(1) == 3 % -> face 134
            nodes_to_remove = sorted_nodes(2); % remove local node 2 (check if local node 2 is bnd or int node)
            if any(sorted_nodes(I_bnd) == nodes_to_remove)
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove);
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove);
            end
        elseif iq_tri2(1) == 4 % -> face 234
            nodes_to_remove = sorted_nodes(1); % remove local node 1 (check if local node 1 is bnd or int node)
            if any(sorted_nodes(I_bnd) == nodes_to_remove)
                bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                bnd_nodes_to_remove = unique(bnd_nodes_to_remove);
            else
                int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                int_nodes_to_remove = unique(int_nodes_to_remove);
            end
        end
    elseif size(I_fix,2) == 1 % there is 1 fixed node
        if I_fix == 1 % local node 1 is fixed -> triangles 123, 124 and 134
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 2 % local node 2 is fixed -> triangles 123, 124 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 3 % local node 3 is fixed -> triangles 123, 134 and 234
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
                % (check if local node 4 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        elseif I_fix == 4 % local node 4 is fixed -> triangles 124, 134 and 234
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
                % (check if local node 3 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
                % (check if local node 2 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            elseif iq_tri2(1) == 3 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
                % (check if local node 1 is bnd or int node)
                if any(sorted_nodes(I_bnd) == nodes_to_remove)
                    bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
                    bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
                else
                    int_nodes_to_remove = [int_nodes_to_remove; nodes_to_remove];
                    int_nodes_to_remove = unique(int_nodes_to_remove(:));
                end
            end
        end
    elseif size(I_fix,2) == 2 % there are 2 fixed nodes
        if sum(ismember([1 2],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            end
        elseif sum(ismember([1 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            end
        elseif sum(ismember([1 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 2 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            end
        elseif sum(ismember([2 3],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 123
                nodes_to_remove = sorted_nodes(4); % remove local node 4
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif sum(ismember([2 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 124
                nodes_to_remove = sorted_nodes(3); % remove local node 3
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        elseif sum(ismember([3 4],I_fix)) == 2
            if iq_tri2(1) == 1 % -> face 134
                nodes_to_remove = sorted_nodes(2); % remove local node 2
            elseif iq_tri2(1) == 2 % -> face 234
                nodes_to_remove = sorted_nodes(1); % remove local node 1
            end
        end
        bnd_nodes_to_remove = [bnd_nodes_to_remove; nodes_to_remove];
        bnd_nodes_to_remove = unique(bnd_nodes_to_remove(:));
    end
end

end % END OF SUBFUNCTION compare_q_slivers_2bnd_nodes

function MESH = choose_case(SETTINGS)
% Usage: MESH = choose_case(SETTINGS)
%
% Purpose:
%   Debug the code by creating a very simple mesh with slivers for
%   different configurations, e.g., nodes in bnd1, bnd2, int, share nodes,
%   share bars, etc.
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%
%   MESH     : [structure] : structure containing the mesh
%
% JMT Jul 2016

prompt = 'Which case? ';
option = input(prompt,'s');
pfix   = [];
pbnd1  = [];
pbnd2  = [];
switch option
    case 'case1'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'             Case 1: Interior sliver with no shared edges\n');
        fprintf(1,'=========================================================================\n');
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 2/3 0.3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case2'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'             Case 2: Interior slivers with 1 shared node\n');
        fprintf(1,'=========================================================================\n');
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 0.7 1/3 0.55; 2/3 2/3 0.45; 1/3 2/3 0.5];
    case 'case3'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'             Case 3: Interior slivers with 1 shared edge\n');
        fprintf(1,'=========================================================================\n');
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 2/3 1/3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case4'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'       Case 4: Boundary sliver (1 bnd node) with no shared edges\n');
        fprintf(1,'=========================================================================\n');
        pbnd2       = [0 0 1];
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 1 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 2/3 0.3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case4a'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'         Case 4a: Boundary sliver (1 bnd node) with 1 shared edge\n');
        fprintf(1,'                 (shared edge contains the bnd node)\n');
        fprintf(1,'=========================================================================\n');
        pbnd2       = [0 0 1];
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 1 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 2/3 1/3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case4b'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'         Case 4b: Boundary sliver (1 bnd node) with 1 shared edge\n');
        fprintf(1,'                 (shared edge does not contain the bnd node)\n');
        fprintf(1,'=========================================================================\n');
        pbnd2       = [1 0 1];
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 1 1; 0 1 1; 1/3 1/3 0.55; 2/3 1/3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case5'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'        Case 5: Boundary sliver (2 bnd nodes) with no shared edges\n');
        fprintf(1,'=========================================================================\n');
        pbnd2       = [0 0 1; 0 1 1];
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 1 0 1; 1 1 1; 1/3 1/3 0.55; 2/3 0.3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case6'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'        Case 6: Boundary sliver (3 bnd nodes) with no shared edges\n');
        fprintf(1,'=========================================================================\n');
        pbnd2       = [0 0 1; 0 1 1; 1/3 1/3 0.55];
        pint        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 1 0 1; 1 1 1; 2/3 0.3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
    case 'case7'
        fprintf(1,'\n\n');
        fprintf(1,'=========================================================================\n');
        fprintf(1,'                     Case 7: nodes in 2 boundaries\n');
        fprintf(1,'=========================================================================\n');
        pbnd1       = [0 0 1; 0 1 1];
        pbnd2       = [1 0 0; 1 1 0];
        pint        = [0 0 0; 0 1 0; 1 0 1; 1 1 1; 1/3 1/3 0.55; 2/3 0.3 0.5; 2/3 2/3 0.5; 1/3 2/3 0.5];
end
MESH.GCOORD = [pfix;pbnd1;pbnd2;pint];
MESH.EL2NOD = delaunay(MESH.GCOORD);
MESH.q      = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
MESH.pfix   = pfix;
MESH.pbnd1  = pbnd1;
MESH.pbnd2  = pbnd2;
MESH.pint   = pint;
if SETTINGS.show_figs
    figure(44)
    clf
    hold on
    axis equal
    view(142.5,30)
    xlabel('X');ylabel('Y');zlabel('Z')
    scatter3(MESH.GCOORD(MESH.EL2NOD,1),MESH.GCOORD(MESH.EL2NOD,2),MESH.GCOORD(MESH.EL2NOD,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
    faceColor = [0.6875 0.8750 0.8984];
    tetramesh(MESH.EL2NOD,MESH.GCOORD,'FaceColor',faceColor,'FaceAlpha',0);
end
end % END OF SUBFUNCTION choose_case

function plot_sphere(r)
% Usage: plot_sphere(r)
%
% Purpose:
%   Plot a sphere centered in the origin of coordinates (0,0,0)
%
% Input:
%   r : [scalar] : radius of the sphere 
%
% Output:
%   none (plot figure)
%
% JMT Jun 2016

lightGrey                    = 0.85*[1 1 1];
[x_sphere,y_sphere,z_sphere] = sphere(50);
x_sphere                     = r*x_sphere;
y_sphere                     = r*y_sphere;
z_sphere                     = r*z_sphere;
surface(x_sphere,y_sphere,z_sphere,'FaceColor', 'none','EdgeColor',lightGrey)
axis([-r r -r r -r r])
axis equal
view(142.5,30)
xlabel('X');ylabel('Y');zlabel('Z')
hold on

end % END OF SUBFUNCTION plot_sphere