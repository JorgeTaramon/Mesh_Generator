function MESH = remove_bad_elements(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: MESH = remove_bad_elements(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   Find bad elements (q < q_bad) and improve them.
%   For every bad element compute the state of stretch/compression for each
%   bar (= relative bar-length change which can be stretch (> 0) or
%   compression (< 0)), and the distortion for each bar 
%   (= abs(relative bar-length change)). 
%   Then do, for each bad element, either:
%   a) Adding a new node in the middle of the longest bar (stretch)
%   b) Rejecting one of the nodes of the shortest bar (compression)
%   This routine takes into account which kind of node (bnd1,bnd2 or int)
%   is being modified (either add or reject). Therefore each element has
%   to be dealt individually with a "for" loop. However, the number of
%   elements with q < q_bad is usually small in comparison with the total
%   number of elements, therefore this routine is not very time consuming.
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
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% DEBUG MODE
%==========================================================================
if nargin == 0
    SETTINGS.q_bad       = 0.25;
    SETTINGS.r_int       = 3471;
    SETTINGS.r_ext       = 6371;
    SETTINGS.h0          = 1;
    SETTINGS.first_guess = 'new';
    SETTINGS.refinement  = 'regular';
    SETTINGS.show_figs   = 1;
    GUIDE_MESH           = struct([]);
    INTERFACE            = struct([]);
    MESH.pfix            = [];
    MESH.pbnd1           = [];
    MESH.pbnd2           = [0  0  0; 1  0  0; 1 1  0;  0 1  0;  0  0 1; 1  0 1; 1 1 1;  0 1 1];
%     MESH.pint            = [1/3 1/3 0.55; 2/3 0.4 0.5; 2/3 2/3 0.5; 1/3 2/3 0.55];
    MESH.pint            = [2/3 1/3 0.55; 2/3 2/3 0.5; 1/3 2/3 0.5; 1/3 1/3 0.55];
%     MESH.pint            = [0.2 0.2 0.2; 0.8 0.8 0.8];
    MESH.GCOORD          = [MESH.pfix;MESH.pbnd1;MESH.pbnd2;MESH.pint];
    MESH.EL2NOD          = delaunay(MESH.GCOORD);
    [MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
    MESH.q               = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
    MESH.rel_change      = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
    figure(212)
        clf
        hold on
        axis equal
        view(155,45)
        xlabel('X');ylabel('Y');zlabel('Z')
        faceColor = [0.6875 0.8750 0.8984];
        tetramesh(MESH.EL2NOD,MESH.GCOORD,'FaceColor',faceColor,'FaceAlpha',0);
        axis([0 1 0 1 0 1])
else
    if SETTINGS.show_figs
        figure(212)
        clf
        hold on
        axis equal
        view(142.5,30)
        xlabel('X');ylabel('Y');zlabel('Z')
        faceColor           = [0.6875 0.8750 0.8984];
        lightGrey           = 0.85*[1 1 1];
        [x_sph,y_sph,z_sph] = sphere(30);
        x_sph               = x_sph*6371;
        y_sph               = y_sph*6371;
        z_sph               = z_sph*6371;
        surface(x_sph,y_sph,z_sph,'FaceColor','none','EdgeColor',lightGrey)
        axis([-6371 6371 -6371 6371 -6371 6371])
        grid on
    end
end

%==========================================================================
% LOAD VARIABLES AND COMPUTE NEEDED INFORMATION
%==========================================================================
q_bad         = SETTINGS.q_bad;
GCOORD        = MESH.GCOORD;
EL2NOD        = MESH.EL2NOD;
pfix          = MESH.pfix;
pbnd1         = MESH.pbnd1;
pbnd2         = MESH.pbnd2;
pint          = MESH.pint;
nfix          = size(MESH.pfix,1);
nbnd1         = size(MESH.pbnd1,1);
nbnd2         = size(MESH.pbnd2,1);
q             = MESH.q;
bars          = MESH.bars;
rel_change    = MESH.rel_change;
barmid_new    = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2;
[q_sorted,iq] = sort(q,'ascend');                  % sort the quality factor of the elements from the worst one to the best one
EL2NOD_sorted = EL2NOD(iq,:);                      % sort the elements from the worst one to the best one
bad_el        = EL2NOD_sorted(q_sorted < q_bad,:); % find sorted bad elements starting from the worst one

%==========================================================================
% COMPUTE NODES TO BE ADDED AND NODES TO BE REJECTED (BND1, BND2, INT)
%==========================================================================
pbnd1_new           = []; % empty pbnd1_new just to enter into the loop
pbnd2_new           = []; % empty pbnd2_new just to enter into the loop
pint_new            = []; % empty pint_new just to enter into the loop
bnd1_node_to_reject = []; % empty bnd1_node_to_reject just to enter into the loop
bnd2_node_to_reject = []; % empty bnd2_node_to_reject just to enter into the loop
int_node_to_reject  = []; % empty int_node_to_reject just to enter into the loop
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    ifix1           = find(abs(sqrt(sum(pfix.^2,2))-SETTINGS.r_int) < 1e-8); % indices for fixed nodes on boundary 1
    ibnd1           = (nfix+1:nfix+nbnd1)';                                  % indices for bnd nodes on boundary 1
else
    ifix1           = [];
    ibnd1           = [];
end
ifix2               = find(abs(sqrt(sum(pfix.^2,2))-SETTINGS.r_ext) < 1e-8); % indices for fixed nodes on boundary 2
ibnd2               = (nfix+nbnd1+1:nfix+nbnd1+nbnd2)';                      % indices for bnd nodes on boundary 2
for i = 1:size(bad_el,1)
    % sort the ith bad element nodes from smallest to biggest (this is done in order to find the bars easily)
    sorted_bad_el               = sort(bad_el(i,:),2); 
    % bad element bars sorted in two columns (each colum is the node at each end of the bar)
    bad_el_Nod2bars             = transpose(reshape(transpose(sorted_bad_el(1,[1 2 1 3 1 4 2 3 2 4 3 4])),2,[])); 
    [~,index]                   = ismember(bars,bad_el_Nod2bars,'rows');
    bad_el_bars                 = find(index); % it might be less than 6. See bar_length.m (bars linking fixed nodes are not included in the bars list since they cannot be modified)
    if isempty(bad_el_bars)
        continue
    end
    bad_el_Nod2bars             = bars(bad_el_bars,:); % take the bars that actually can be modified
    rel_change_bad_el           = rel_change(bad_el_bars);
    rel_change_distorted_bad_el = abs(rel_change_bad_el);
    [~,I]                       = sort(rel_change_distorted_bad_el,'descend');
    % Take the most distorted bar of each bad element and check if it is stretched (> 0) or compressed (< 0)
        % ADD NODE
        if rel_change_bad_el(I(1)) > 0
            % We have to distinguish between added nodes on boundary 1, boundary 2 or in the interior
            % BND1: the ends of the bar are on the bnd1
            if (sum(bad_el_Nod2bars(I(1),1) == ifix1) || sum(bad_el_Nod2bars(I(1),1) == ibnd1)) && ...
               (sum(bad_el_Nod2bars(I(1),2) == ifix1) || sum(bad_el_Nod2bars(I(1),2) == ibnd1))
                pbnd1_new = [pbnd1_new; barmid_new(bad_el_bars(I(1)),:)];
            % BND2: the ends of the bar are on the bnd2
            elseif (sum(bad_el_Nod2bars(I(1),1) == ifix2) || sum(bad_el_Nod2bars(I(1),1) == ibnd2)) && ...
                   (sum(bad_el_Nod2bars(I(1),2) == ifix2) || sum(bad_el_Nod2bars(I(1),2) == ibnd2))
                pbnd2_new = [pbnd2_new; barmid_new(bad_el_bars(I(1)),:)];
            % INT: at least one of the ends of the bar is in the interior
            else
                pint_new = [pint_new; barmid_new(bad_el_bars(I(1)),:)]; 
            end
            if SETTINGS.show_figs
                figure(212)
                hold on
                grid on
%                 scatter3(GCOORD(sorted_bad_el(1,:),1),GCOORD(sorted_bad_el(1,:),2),GCOORD(sorted_bad_el(1,:),3), ...
%                     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
                tetramesh(sorted_bad_el,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
                if nargin == 0
                    figure(i)
                    tetramesh(sorted_bad_el,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
                    axis equal
                    view(155,45)
                    xlabel('X');ylabel('Y');zlabel('Z')
                    axis([0 1 0 1 0 1])
                    grid on
                end
            end
        % REJECT NODE (we have to specify between boundary nodes and interior nodes)
        else
            % reject bnd1 node
            if (bars(bad_el_bars(I(1)),2) > nfix && bars(bad_el_bars(I(1)),2) <= nfix+nbnd1)
                bnd1_node_to_reject = [bnd1_node_to_reject; bars(bad_el_bars(I(1)),2)];
            % reject bnd2 node
            elseif bars(bad_el_bars(I(1)),2) > nfix+nbnd1 && bars(bad_el_bars(I(1)),2) <= nfix+nbnd1+nbnd2
                bnd2_node_to_reject = [bnd2_node_to_reject; bars(bad_el_bars(I(1)),2)];
            % reject int node
            elseif bars(bad_el_bars(I(1)),2) > nfix+nbnd1+nbnd2
                int_node_to_reject = [int_node_to_reject; bars(bad_el_bars(I(1)),2)];
            end
            if SETTINGS.show_figs
                figure(212)
                hold on
                grid on
%                 scatter3(GCOORD(sorted_bad_el(1,:),1),GCOORD(sorted_bad_el(1,:),2),GCOORD(sorted_bad_el(1,:),3), ...
%                     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
                tetramesh(sorted_bad_el,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
                if nargin == 0
                    figure(i)
                    tetramesh(sorted_bad_el,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
                    axis equal
                    view(155,45)
                    xlabel('X');ylabel('Y');zlabel('Z')
                    axis([0 1 0 1 0 1])
                    grid on
                end
            end
        end
end

%===============================================================================
% REMOVE REPEATED NEW NODES AND PROJECT THE ADDED BND NODES ONTO THE BOUNDARIES 
%===============================================================================
% BOUNDARY 1
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    pbnd1_new = unique(pbnd1_new,'rows','stable');   % remove repeated pbnd1_new
    if ~isempty(pbnd1_new)
        % Project these new bnd1 nodes onto the spherical boundary1
        pbnd1_new_SPH      = cartesian2spherical(pbnd1_new);               % convert to spherical coordinates
        pbnd1_new_SPH(:,3) = SETTINGS.r_int*ones(size(pbnd1_new_SPH,1),1); % change the r-coordinate by r_int (projecting onto bnd 1)
        pbnd1_new          = spherical2cartesian(pbnd1_new_SPH);           % Cartesian coordinates of new bnd1 nodes after projecting on bnd1
    end
else
    pbnd1_new = [];
end
% BOUNDARY 2
pbnd2_new = unique(pbnd2_new,'rows','stable'); % remove repeated pbnd2_new
if ~isempty(pbnd2_new)
    % Project these new bnd2 nodes onto the spherical boundary2
    pbnd2_new_SPH      = cartesian2spherical(pbnd2_new);               % convert to spherical coordinates
    pbnd2_new_SPH(:,3) = SETTINGS.r_ext*ones(size(pbnd2_new_SPH,1),1); % change the r-coordinate by r_ext (projecting onto bnd 2)
    pbnd2_new          = spherical2cartesian(pbnd2_new_SPH);           % Cartesian coordinates of new bnd2 nodes after projecting on bnd2
end
% INTERIOR
pint_new = unique(pint_new,'rows','stable'); % remove repeated pbnd1_new

%==========================================================================
% REMOVE REPEATED NODES TO REJECT
%==========================================================================
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    bnd1_node_to_reject = unique(bnd1_node_to_reject,'stable'); % remove repeated bnd1_node_to_reject
else
    bnd1_node_to_reject = [];
end
bnd2_node_to_reject = unique(bnd2_node_to_reject,'stable'); % remove repeated bnd2_node_to_reject
int_node_to_reject  = unique(int_node_to_reject,'stable');  % remove repeated int_node_to_reject

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(41)
    clf
    hold on
    grid on
    if ~isempty(pbnd1_new)
        scatter3(pbnd1_new(:,1),pbnd1_new(:,2),pbnd1_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd1 node added in yellow
    end
    if ~isempty(pbnd2_new)
        scatter3(pbnd2_new(:,1),pbnd2_new(:,2),pbnd2_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd2 node added in yellow
    end
    if ~isempty(pint_new)
        scatter3(pint_new(:,1),pint_new(:,2),pint_new(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]) % int node added in green
    end
    if ~isempty(bnd1_node_to_reject)
        scatter3(GCOORD(bnd1_node_to_reject,1),GCOORD(bnd1_node_to_reject,2),GCOORD(bnd1_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) % bnd1 node rejected in blue
    end
    if ~isempty(bnd2_node_to_reject)
        scatter3(GCOORD(bnd2_node_to_reject,1),GCOORD(bnd2_node_to_reject,2),GCOORD(bnd2_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) % bnd2 node rejected in blue
    end
    if ~isempty(int_node_to_reject)
        scatter3(GCOORD(int_node_to_reject,1),GCOORD(int_node_to_reject,2),GCOORD(int_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]) % int node rejected in red
    end
    axis equal
    view(142.5,30)
    xlabel('X');ylabel('Y');zlabel('Z')
    title({'Bad elements in blue';'bnd nodes: added (yellow), rejected (blue); int nodes: added (green), rejected (red)'})
end

%==========================================================================
% REJECT AND ADD THE NODES COMPUTED BEFORE
%==========================================================================
if SETTINGS.r_int > 0 && ~isempty(pbnd1)
    pbnd1(bnd1_node_to_reject-nfix,:) = []; % reject the 2nd node of the most compressed bar of the bad element (bnd1 nodes)
    pbnd1                             = [pbnd1; pbnd1_new]; % add new bnd1 nodes
end
pbnd2(bnd2_node_to_reject-nfix-nbnd1,:)     = []; % reject the 2nd node of the most compressed bar of the bad element (bnd2 nodes)
pbnd2                                       = [pbnd2; pbnd2_new]; % add new bnd2 nodes
pint(int_node_to_reject-nfix-nbnd1-nbnd2,:) = []; % reject the 2nd node of the most compressed bar of the bad element (int nodes)
pint                                        = [pint; pint_new];   % add new int nodes
GCOORD_new                                  = [pfix; pbnd1; pbnd2; pint];
EL2NOD_new                                  = delaunay(GCOORD_new);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH  = cartesian2spherical(GCOORD_new);
EL2NOD_new  = EL2NOD_new(~(sum(ismember(EL2NOD_new,find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
q_new       = tetra_mesh_quality(GCOORD_new,EL2NOD_new);
new_bad_el  = EL2NOD_new(q_new < q_bad,:);

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd1                     = pbnd1;
MESH.pbnd2                     = pbnd2;
MESH.pint                      = pint;
MESH.GCOORD                    = GCOORD_new;
MESH.EL2NOD                    = EL2NOD_new;
MESH.q                         = q_new;
MESH.s                         = shape_measure(MESH.GCOORD,MESH.EL2NOD);
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
MESH.rel_change                = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
MESH.rel_change_abs            = abs(MESH.rel_change);        % absolute value of relative bar-length change
MESH.mean_misfit_bar_length    = sum(MESH.rel_change_abs)/size(MESH.L,1);

%==========================================================================
% DEBUG MODE
%==========================================================================
if nargin == 0
    if SETTINGS.show_figs
        figure(410)
        clf
        scatter3(GCOORD_new(EL2NOD_new,1),GCOORD_new(EL2NOD_new,2),GCOORD_new(EL2NOD_new,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
        hold on
        tetramesh(EL2NOD_new,GCOORD_new,'FaceColor',faceColor,'FaceAlpha',0);
        tetramesh(new_bad_el,GCOORD_new,'FaceColor',faceColor,'FaceAlpha',0.3);
        if ~isempty(pbnd1_new)
            scatter3(pbnd1_new(:,1),pbnd1_new(:,2),pbnd1_new(:,3), ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd1 node added in yellow
        elseif ~isempty(pbnd2_new)
            scatter3(pbnd2_new(:,1),pbnd2_new(:,2),pbnd2_new(:,3), ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % bnd1 node added in yellow
        elseif ~isempty(pint_new)
            scatter3(pint_new(:,1),pint_new(:,2),pint_new(:,3), ...
                'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]) % bnd1 node added in green
        end
        scatter3(GCOORD(bnd1_node_to_reject,1),GCOORD(bnd1_node_to_reject,2),GCOORD(bnd1_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1]) % bnd1 node rejected in blue
        scatter3(GCOORD(bnd2_node_to_reject,1),GCOORD(bnd2_node_to_reject,2),GCOORD(bnd2_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1]) % bnd2 node rejected in blue
        scatter3(GCOORD(int_node_to_reject,1),GCOORD(int_node_to_reject,2),GCOORD(int_node_to_reject,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]) % int node rejected in red
        axis equal
        view(155,45)
        xlabel('X');ylabel('Y');zlabel('Z')
        title('mesh after remove bad elements')
    end
end
  
end % END OF FUNCTION remove_bad_elements