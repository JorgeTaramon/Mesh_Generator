function [pfix,pbnd,segBC,pseg] = boundary_nodes(SETTINGS,GUIDE_MESH)
% Usage: [pfix,pbnd,segBC,pseg] = boundary_nodes(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Generation of boundary nodes for a 2D rectangular mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   pfix       : [matrix]    : coord of fixed nodes
%   pbnd       : [matrix]    : coord of boundary nodes
%   segBC      : [matrix]    : info of boundary segments
%   pseg       : [vector]    : pointID for bnd nodes
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth      = SETTINGS.depth;         % domain depth (km)
length     = SETTINGS.length;        % domain length (km)
x0         = GUIDE_MESH.x0;          % point around which the refined and transition zones are defined
z0         = GUIDE_MESH.z0;          % point around which the refined and transition zones are defined
l_tran     = GUIDE_MESH.l_tran;      % length of transition zone (km)
x_tran_l   = x0 - l_tran/2;          % left boundary of the transition zone
x_tran_r   = x0 + l_tran/2;          % right boundary of the transition zone
l_ref      = GUIDE_MESH.l_ref;       % length of refined zone (km)
x_ref_l    = x0 - l_ref/2;           % left boundary of the refined zone
x_ref_r    = x0 + l_ref/2;           % right boundary of the refined zone

%==========================================================================
% DEFINE PFIX
%==========================================================================
pfix = [x0 - length/2  z0 - depth ; ...
        x0 + length/2  z0 - depth ; ...
        x0 + length/2  z0 ; ...
        x0 - length/2  z0 ]; % vertex points COUNTER-CLOCKWISE
nfix = size(pfix,1);
nseg = nfix;

%==========================================================================
% DEFINE EACH BOUNDARY SEGMENT
%==========================================================================
seg_connect    = [(1:nseg)',[(2:nseg)';1]];  % This is the connectivity matrix for the fixed points
segBC          = zeros(nseg,6);              % pre-allocate memory
segBC(:,[1 2]) = pfix(seg_connect(:,1),:);   % p0 = (x0,y0) for seg_start
segBC(:,[3 4]) = pfix(seg_connect(:,2),:);   % p1 = (x1,y1) for seg_end
segBC(:,[5 6]) = GetLineSlope(segBC(:,[1 2]),segBC(:,[3 4]),nseg); % slope_info = [a,b]

%==========================================================================
% COMPUTE PBND
%==========================================================================
pbnd = [];
pseg = [];
for i = 1:nseg

    if i == 1
        j = 2;
        k = 2;
    elseif i == 2
        j = 1;
        k = 1;
    elseif i == 3
        j = 2;
        k = 4;
    elseif i == 4
        j = 1;
        k = 3;
    end
    
    % COARSE ZONE
    pbnd_this_segment_coarse = regular_bnd_nodes(segBC,GUIDE_MESH.l0_coarse,i);
    x_this_segment_coarse    = pbnd_this_segment_coarse(:,1);
    
    if sum(ismember([GUIDE_MESH.l0_ref; GUIDE_MESH.l0_coarse],unique(GUIDE_MESH.L0_guide(GUIDE_MESH.GCOORD_GUIDE(:,j) == segBC(i,k))))) == 2
        
        pbnd_this_segment_coarse(x_this_segment_coarse > x_tran_l &...
                                 x_this_segment_coarse < x_tran_r,:) = []; % remove bnd nodes that are in the transition zone
        pbnd_this_segment_coarse = reject_points_line(pbnd_this_segment_coarse,GUIDE_MESH); % reject pbnd_this_segment_coarse using probability
        
        % TRANSITION ZONE
        pbnd_this_segment_tran   = regular_bnd_nodes(segBC,GUIDE_MESH.l0_ref,i);
        x_this_segment_tran      = pbnd_this_segment_tran(:,1);
        pbnd_this_segment_tran   = pbnd_this_segment_tran(x_this_segment_tran > x_tran_l &...
                                                          x_this_segment_tran < x_tran_r,:); % take only the nodes in the transition zone
        x_tran_temp              = x_this_segment_tran(x_this_segment_tran > x_tran_l & ...
                                                       x_this_segment_tran < x_tran_r);  % take x inside the transition zone
        pbnd_this_segment_tran(x_tran_temp > x_ref_l & ...
                               x_tran_temp < x_ref_r,:) = []; % remove bnd nodes that are inside the refined zone
        pbnd_this_segment_tran   = reject_points_line(pbnd_this_segment_tran,GUIDE_MESH); % reject pbnd_this_segment_tran using probability

        % REFINED ZONE
        pbnd_this_segment_ref    = regular_bnd_nodes(segBC,GUIDE_MESH.l0_ref,i);
        x_this_segment_ref       = pbnd_this_segment_ref(:,1);
        pbnd_this_segment_ref    = pbnd_this_segment_ref(x_this_segment_ref > x_ref_l &...
                                                         x_this_segment_ref < x_ref_r,:); % take only the nodes in the refined zone
        pbnd_this_segment_ref    = reject_points_line(pbnd_this_segment_ref,GUIDE_MESH); % reject pbnd_this_segment_ref using probability
    else
        pbnd_this_segment_coarse = reject_points_line(pbnd_this_segment_coarse,GUIDE_MESH); % reject pbnd_this_segment_coarse using probability
        pbnd_this_segment_tran   = [];
        pbnd_this_segment_ref    = [];
    end
    pbnd = [pbnd; pbnd_this_segment_coarse; pbnd_this_segment_tran; pbnd_this_segment_ref];
    pseg = [pseg i*ones(1,size([pbnd_this_segment_coarse; pbnd_this_segment_tran; pbnd_this_segment_ref],1))];
end

end % END OF FUNCTION boundary_nodes

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function GCOORD = reject_points_line(GCOORD,GUIDE_MESH)
% Usage: GCOORD = reject_points_line(GCOORD,GUIDE_MESH)
%
% Purpose:
%   Reject points on a line using point density (based on Persson and
%   Strang, 2004)
%
% Input:
%   GCOORD     : [matrix]    : coordinates of points
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   GCOORD     : [matrix]    : coordinates of points
%
% JMT Jun 2017

if isempty(GCOORD)
    return
end

L0_GCOORD = bar_L0_guide(GCOORD,GUIDE_MESH); % compute the desired length for each point
r0_GCOORD = sqrt(2)./L0_GCOORD;              % probability to keep the points (it is proportional to L^-1 because points are on a line)
GCOORD    = GCOORD(rand(size(GCOORD,1),1)<r0_GCOORD./max(r0_GCOORD),:); % reject points with a probability proportional to L^-1

end % END OF FUNCTION reject_points_line