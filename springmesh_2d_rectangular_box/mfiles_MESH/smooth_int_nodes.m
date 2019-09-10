function [MESH] = smooth_int_nodes(MESH,GUIDE_MESH,SETTINGS)
% Usage: [MESH] = smooth_int_nodes(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose: 
%   For each element which quality factor is smaller than a quality
%   tolerance (q < q_smooth), select the interior nodes that belong to
%   those elements. Then, compute the average of the position of the
%   centres of the elements that surround each interior node. Finally, move
%   each interior node to the position given by the average of the centres
%   of elements that surround it
%
% Input:
%   MESH       : [structure] : structure containing the mesh
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   MESH       : [structure] : structure containing the mesh
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
q_smooth = SETTINGS.q_smooth;
EL2NOD   = MESH.EL2NOD;
GCOORD   = MESH.GCOORD;
pint     = MESH.pint;
nfix     = size(MESH.pfix,1);
nbnd     = size(MESH.pbnd,1);

%==========================================================================
% COMPUTE NEW POSITIONS FOR PINT SMOOTHED
%==========================================================================
[q_sorted,I]             = sort(MESH.q); % sort q from the worst q to the best q
EL2NOD_sorted            = EL2NOD(I(q_sorted < q_smooth),:); % select those elements with q < q_smooth from the worst one to the best one
int_nodes_to_be_smoothed = unique(EL2NOD_sorted(EL2NOD_sorted > nfix+nbnd),'stable'); % select the interior nodes that will be smoothed
if size(int_nodes_to_be_smoothed,1) == 1
    int_nodes_to_be_smoothed = int_nodes_to_be_smoothed';
end
% preallocated vectors (for x and y coordinates) for the average of all element centers that surround each interior node
x_average_all_centres_surround_each_int_node = zeros(size(int_nodes_to_be_smoothed,1),1); 
y_average_all_centres_surround_each_int_node = zeros(size(int_nodes_to_be_smoothed,1),1);

x_centre_plot = [];
y_centre_plot = [];

for i = 1:size(int_nodes_to_be_smoothed,1) % counting the interior nodes
    % boolean indexes to know which element contains the i-th interior node
    b                                               = reshape(EL2NOD(:) == int_nodes_to_be_smoothed(i),size(EL2NOD,1),3);
    % conectivity matrix with the elements that contain the i-th interior node
    d                                               = EL2NOD(sum(b,2)==1,:);
    % x coordinate for the average of all element centers that surround the i-th interior node
    x_average_el_centre_surround_this_int_node      = mean(reshape(GCOORD(d(:),1),size(d,1),3),2);
    % y coordinate for the average of all element centers that surround the i-th interior node
    y_average_el_centre_surround_this_int_node      = mean(reshape(GCOORD(d(:),2),size(d,1),3),2);
    % x coordinate for the average of all element centers that surround each interior node
    x_average_all_centres_surround_each_int_node(i) = mean(x_average_el_centre_surround_this_int_node);
    % y coordinate for the average of all element centers that surround each interior node    
    y_average_all_centres_surround_each_int_node(i) = mean(y_average_el_centre_surround_this_int_node);
    x_centre_plot                                   = [x_centre_plot; x_average_el_centre_surround_this_int_node];
    y_centre_plot                                   = [y_centre_plot; y_average_el_centre_surround_this_int_node];
end

pint_to_be_modified  = [x_average_all_centres_surround_each_int_node ...
                        y_average_all_centres_surround_each_int_node];

pint_smoothed                                       = pint;
pint_smoothed(int_nodes_to_be_smoothed-nfix-nbnd,:) = pint_to_be_modified;

%==========================================================================
% DATA FOR PLOTS AND OUTPUT
%==========================================================================
GCOORD_smoothed                            = [MESH.pfix;MESH.pbnd;pint_smoothed];
EL2NOD_smoothed                            = delaunay(GCOORD_smoothed);
q_smoothed                                 = mesh_quality(GCOORD_smoothed,EL2NOD_smoothed);
MESH_SMOOTHED.GCOORD                       = GCOORD_smoothed;
MESH_SMOOTHED.EL2NOD                       = EL2NOD_smoothed;
MESH_SMOOTHED.pfix                         = MESH.pfix;
[L0_smoothed,L_smoothed,bars_smoothed,~,~] = bar_length(MESH_SMOOTHED,GUIDE_MESH,SETTINGS);
MESH_SMOOTHED.rel_change                   = (L_smoothed - L0_smoothed)./L0_smoothed; % relative bar-length change
MESH_SMOOTHED.rel_change_abs               = abs(MESH_SMOOTHED.rel_change); % absolute value of relative bar-length change
MESH_SMOOTHED.mean_misfit_bar_length       = sum(MESH_SMOOTHED.rel_change_abs)/size(L_smoothed,1);
MESH_SMOOTHED.iter                         = MESH.iter;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(8)
    clf
    subplot(1,2,1)
    trimesh(EL2NOD_sorted,GCOORD(:,1),GCOORD(:,2),'Color',[0 0 0]) % plot mesh before smoothing
    hold on
    line([GCOORD(int_nodes_to_be_smoothed,1) pint_to_be_modified(:,1)]',...
         [GCOORD(int_nodes_to_be_smoothed,2) pint_to_be_modified(:,2)]','Color',[0 1 0])
    plot(GCOORD(int_nodes_to_be_smoothed,1),GCOORD(int_nodes_to_be_smoothed,2),'rx')
    plot(pint_to_be_modified(:,1),pint_to_be_modified(:,2),'gx')
    axis equal
    title('before smoothing')
    
    subplot(1,2,2)
    trimesh(EL2NOD_sorted,GCOORD_smoothed(:,1),GCOORD_smoothed(:,2),'Color',[0 0 0]) % plot mesh before smoothing
    hold on
    line([GCOORD(int_nodes_to_be_smoothed,1) pint_to_be_modified(:,1)]',...
         [GCOORD(int_nodes_to_be_smoothed,2) pint_to_be_modified(:,2)]','Color',[0 1 0])
    plot(GCOORD(int_nodes_to_be_smoothed,1),GCOORD(int_nodes_to_be_smoothed,2),'rx')
    plot(pint_to_be_modified(:,1),pint_to_be_modified(:,2),'gx')
    axis equal
    title('after smoothing')
    
    figure(9)
    clf
    subplot(2,2,1)
    tmp                = SETTINGS.save_figs;
    SETTINGS.save_figs = 0; % don't save the subplots
    plot_q_factor_histogram(MESH.q,MESH.iter,SETTINGS)
    title('Histogram of quality factor before smoothing')
        
    subplot(2,2,2)
    plot_q_factor_histogram(q_smoothed,MESH.iter,SETTINGS)
    title('Histogram of quality factor after smoothing')
    
    subplot(2,2,3)
    plot_misfit_bar_length(MESH,SETTINGS)
    title('Histogram of misfit bar-length before smoothing')
    
    subplot(2,2,4)
    plot_misfit_bar_length(MESH_SMOOTHED,SETTINGS)
    title('Histogram of misfit bar-length after smoothing')
    
    SETTINGS.save_figs = tmp;
end

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.q                      = mesh_quality(GCOORD_smoothed,EL2NOD_smoothed);
MESH.pint                   = pint_smoothed;
MESH.GCOORD                 = GCOORD_smoothed;
MESH.EL2NOD                 = EL2NOD_smoothed;
MESH.L                      = L_smoothed;
MESH.L0                     = L0_smoothed;
MESH.bars                   = bars_smoothed;
MESH.rel_change             = MESH_SMOOTHED.rel_change;
MESH.rel_change_abs         = MESH_SMOOTHED.rel_change_abs;
MESH.mean_misfit_bar_length = MESH_SMOOTHED.mean_misfit_bar_length;

end % END OF FUNCTION smooth_int_nodes