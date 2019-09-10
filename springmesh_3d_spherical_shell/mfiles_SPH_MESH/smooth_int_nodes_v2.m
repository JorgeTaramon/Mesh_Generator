function MESH = smooth_int_nodes_v2(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: MESH = smooth_int_nodes_v2(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose: 
%   For each element which quality factor is smaller than a quality
%   tolerance (q < q_smooth), select the interior nodes that belong to
%   those elements. Then, compute the average of the element centres that
%   surround each interior node. Finally, move each interior node to the
%   position given by the average of the elements centres that surround
%   them.
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
% JMT Jun 2015
% JMT Jun 2016: cleaned up
% JMT Aug 2017: Handles internal boundaries
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% DEBUG MODE
%==========================================================================
if nargin == 0
    SETTINGS.q_smooth                   = 0.2;
    SETTINGS.r_int                      = 3471;
    SETTINGS.h0                         = 1;
    SETTINGS.mean_misfit_bar_length_tol = 0.14;
    SETTINGS.first_guess                = 'new';
    SETTINGS.refinement                 = 'regular';
    SETTINGS.show_figs                  = 1;
    SETTINGS.save_figs                  = 0;
    GUIDE_MESH                          = struct([]);
    INTERFACE                           = struct([]);
    MESH.iter                           = 0;
    MESH.pfix                           = [];
    MESH.pbnd1                          = [];
    MESH.pbnd2                          = [0  0  0; 1  0  0; 1 1  0;  0 1  0;  0  0 1; 1  0 1; 1 1 1;  0 1 1];
    MESH.pint                           = [0.2 0.2 0.2];
%     MESH.pint                            = [0.2 0.2 0.2; 0.8 0.8 0.8];
    MESH.GCOORD                         = [MESH.pfix;MESH.pbnd1;MESH.pbnd2;MESH.pint];
    MESH.EL2NOD                         = delaunay(MESH.GCOORD);
    [MESH.L0,MESH.L,MESH.bars,~,~]      = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
    MESH.q                              = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
    MESH.rel_change                     = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
    MESH.rel_change_abs                 = abs(MESH.rel_change);        % absolute value of relative bar-length change
    MESH.mean_misfit_bar_length         = sum(MESH.rel_change_abs)/size(MESH.L,1);
end

%==========================================================================
% LOAD VARIABLES
%==========================================================================
q_smooth          = SETTINGS.q_smooth;
EL2NOD            = MESH.EL2NOD;
GCOORD            = MESH.GCOORD;
pint              = MESH.pint;
nfix              = size(MESH.pfix,1);
nbnd1             = size(MESH.pbnd1,1);
nbnd2             = size(MESH.pbnd2,1);
nbnd_ref_face_bot = size(MESH.pbnd_ref_face_bot,1);
nbnd_ref_face_1   = size(MESH.pbnd_ref_face_1,1);
nbnd_ref_face_2   = size(MESH.pbnd_ref_face_2,1);
nbnd_ref_face_3   = size(MESH.pbnd_ref_face_3,1);
nbnd_ref_face_4   = size(MESH.pbnd_ref_face_4,1);

%==========================================================================
% COMPUTE NEW POSITIONS FOR PINT SMOOTHED
%==========================================================================
[q_sorted,I]             = sort(MESH.q); % sort q from the worst q to the best q
EL2NOD_sorted            = EL2NOD(I(q_sorted < q_smooth),:); % select those elements with q < q_smooth from the worst one to the best one
int_nodes_to_be_smoothed = unique(EL2NOD_sorted(EL2NOD_sorted > nfix              + ...
                                                                nbnd1             + ...
                                                                nbnd2             + ...
                                                                nbnd_ref_face_bot + ...
                                                                nbnd_ref_face_1   + ...
                                                                nbnd_ref_face_2   + ...
                                                                nbnd_ref_face_3   + ...
                                                                nbnd_ref_face_4),'stable'); % select the interior nodes that will be smoothed
if size(int_nodes_to_be_smoothed,1) == 1
    int_nodes_to_be_smoothed = int_nodes_to_be_smoothed';
end
% preallocated vectors (for x and y coordinates) for the average of all element centers that surround each interior node
x_average_all_centres_surround_each_int_node = zeros(size(int_nodes_to_be_smoothed,1),1); 
y_average_all_centres_surround_each_int_node = zeros(size(int_nodes_to_be_smoothed,1),1);
z_average_all_centres_surround_each_int_node = zeros(size(int_nodes_to_be_smoothed,1),1);

x_centre_plot = [];
y_centre_plot = [];
z_centre_plot = [];

for i = 1:size(int_nodes_to_be_smoothed,1) % counting the interior nodes
    % boolean indexes to know which element contains the i-th interior node
    b                                               = reshape(EL2NOD(:) == int_nodes_to_be_smoothed(i),size(EL2NOD,1),4);
    % conectivity matrix with the elements that contain the i-th interior node
    d                                               = EL2NOD(sum(b,2)==1,:);
    % x,y,z coordinate for the average of all element centers that surround the i-th interior node
    x_average_el_centre_surround_this_int_node      = mean(reshape(GCOORD(d(:),1),size(d,1),4),2);
    y_average_el_centre_surround_this_int_node      = mean(reshape(GCOORD(d(:),2),size(d,1),4),2);
    z_average_el_centre_surround_this_int_node      = mean(reshape(GCOORD(d(:),3),size(d,1),4),2);
    % x,y,z coordinate for the average of all element centers that surround each interior node
    x_average_all_centres_surround_each_int_node(i) = mean(x_average_el_centre_surround_this_int_node);
    y_average_all_centres_surround_each_int_node(i) = mean(y_average_el_centre_surround_this_int_node);
    z_average_all_centres_surround_each_int_node(i) = mean(z_average_el_centre_surround_this_int_node);
    
    x_centre_plot                                   = [x_centre_plot; x_average_el_centre_surround_this_int_node];
    y_centre_plot                                   = [y_centre_plot; y_average_el_centre_surround_this_int_node];
    z_centre_plot                                   = [z_centre_plot; z_average_el_centre_surround_this_int_node];
end

pint_to_be_modified  = [x_average_all_centres_surround_each_int_node ...
                        y_average_all_centres_surround_each_int_node ...
                        z_average_all_centres_surround_each_int_node];

pint_smoothed                                              = pint;
pint_smoothed(int_nodes_to_be_smoothed - ...
              nfix                     - ...
              nbnd1                    - ...
              nbnd2                    - ...
              nbnd_ref_face_bot        - ...
              nbnd_ref_face_1          - ...
              nbnd_ref_face_2          - ...
              nbnd_ref_face_3          - ...
              nbnd_ref_face_4,:) = pint_to_be_modified;

%==========================================================================
% DATA FOR PLOTS
%==========================================================================
GCOORD_smoothed     = [MESH.pfix;             ...
                      MESH.pbnd1;             ...
                      MESH.pbnd2;             ...
                      MESH.pbnd_ref_face_bot; ...
                      MESH.pbnd_ref_face_1;   ...
                      MESH.pbnd_ref_face_2;   ...
                      MESH.pbnd_ref_face_3;   ...
                      MESH.pbnd_ref_face_4;   ...
                      pint_smoothed];
EL2NOD_smoothed     = delaunay(GCOORD_smoothed);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH_smoothed = cartesian2spherical(GCOORD_smoothed);
EL2NOD_smoothed     = EL2NOD_smoothed(~(sum(ismember(EL2NOD_smoothed,find(abs(GCOORD_SPH_smoothed(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
q_smoothed                                 = tetra_mesh_quality(GCOORD_smoothed,EL2NOD_smoothed);
MESH_SMOOTHED.GCOORD                       = GCOORD_smoothed;
MESH_SMOOTHED.EL2NOD                       = EL2NOD_smoothed;
MESH_SMOOTHED.pfix                         = MESH.pfix;
[L0_smoothed,L_smoothed,bars_smoothed,~,~] = bar_length(MESH_SMOOTHED,GUIDE_MESH,INTERFACE,SETTINGS);
MESH_SMOOTHED.rel_change                   = (L_smoothed - L0_smoothed)./L0_smoothed; % relative bar-length change
MESH_SMOOTHED.rel_change_abs               = abs(MESH_SMOOTHED.rel_change); % absolute value of relative bar-length change
MESH_SMOOTHED.mean_misfit_bar_length       = sum(MESH_SMOOTHED.rel_change_abs)/size(L_smoothed,1);
MESH_SMOOTHED.iter                         = MESH.iter;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(42)
    clf
    subplot(1,2,1)
    if nargin ~= 0
        plot_sphere(SETTINGS.r_ext)
    end
    view(142.5,30)
    xlabel('X');ylabel('Y');zlabel('Z')
    hold on
    scatter3(GCOORD(int_nodes_to_be_smoothed,1),GCOORD(int_nodes_to_be_smoothed,2),GCOORD(int_nodes_to_be_smoothed,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    faceColor = [1 0 0];
    tetramesh(EL2NOD_sorted,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.2);
    scatter3(pint_to_be_modified(:,1),pint_to_be_modified(:,2),pint_to_be_modified(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
    line([GCOORD(int_nodes_to_be_smoothed,1) pint_to_be_modified(:,1)]', ...
         [GCOORD(int_nodes_to_be_smoothed,2) pint_to_be_modified(:,2)]', ...
         [GCOORD(int_nodes_to_be_smoothed,3) pint_to_be_modified(:,3)]','Color',[0 1 0])
    title('before smoothing')
    axis equal
    subplot(1,2,2)
    if nargin ~= 0
        plot_sphere(SETTINGS.r_ext)
    end
    view(142.5,30)
    xlabel('X');ylabel('Y');zlabel('Z')
    hold on
    scatter3(GCOORD(int_nodes_to_be_smoothed,1),GCOORD(int_nodes_to_be_smoothed,2),GCOORD(int_nodes_to_be_smoothed,3), ...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    faceColor = [0.6875 0.8750 0.8984];
    tetramesh(EL2NOD_sorted,GCOORD_smoothed,'FaceColor',faceColor,'FaceAlpha',0.3);
    scatter3(pint_to_be_modified(:,1),pint_to_be_modified(:,2),pint_to_be_modified(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
    line([GCOORD(int_nodes_to_be_smoothed,1) pint_to_be_modified(:,1)]', ...
         [GCOORD(int_nodes_to_be_smoothed,2) pint_to_be_modified(:,2)]', ...
         [GCOORD(int_nodes_to_be_smoothed,3) pint_to_be_modified(:,3)]','Color',[0 1 0])
    title('after smoothing')
    axis equal
    
    figure(43)
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
MESH.q                      = q_smoothed;
MESH.s                      = shape_measure(GCOORD_smoothed,EL2NOD_smoothed);
MESH.pint                   = pint_smoothed;
MESH.GCOORD                 = GCOORD_smoothed;
MESH.EL2NOD                 = EL2NOD_smoothed;
MESH.L                      = L_smoothed;
MESH.L0                     = L0_smoothed;
MESH.bars                   = bars_smoothed;
MESH.rel_change             = MESH_SMOOTHED.rel_change;
MESH.rel_change_abs         = MESH_SMOOTHED.rel_change_abs;
MESH.mean_misfit_bar_length = MESH_SMOOTHED.mean_misfit_bar_length;

end % END OF FUNCTION smooth_int_nodes_v2

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

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

end % END OF SUBFUNCTION plot_sphere