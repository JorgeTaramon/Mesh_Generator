function MESH = pfix_refined_region(MESH,GUIDE_MESH)
% Usage: MESH = pfix_refined_region(MESH,GUIDE_MESH)
%
% Purpose:
%   Generation of pfix for embedded high resolution sub-region
%
% Input:
%   MESH       : [structure] : structure containing the mesh
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   MESH       : [structure] : structure containing the mesh
%
% JMT Jul 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
x0         = GUIDE_MESH.x0;          % point around which the refined and transition zones are defined
z0         = GUIDE_MESH.z0;          % point around which the refined and transition zones are defined
l_ref      = GUIDE_MESH.l_ref;       % length of refined zone (km)
d_ref      = GUIDE_MESH.d_ref;       % max depth in the refined zone (km)
x_ref_l    = x0 - l_ref/2;           % left boundary of the refined zone
x_ref_r    = x0 + l_ref/2;           % right boundary of the refined zone
l0_ref     = GUIDE_MESH.l0_ref;      % desired length (km) for refined zone (constant)
pbnd       = MESH.pbnd;
pseg       = MESH.pseg;

%==========================================================================
% REMOVE PBND ON SURFACE INSIDE THE REFINED REGION
%==========================================================================
x = pbnd(:,1);
z = pbnd(:,2);
pbnd(x >= x_ref_l & x <= x_ref_r & z == z0,:) = [];
pseg(x >= x_ref_l & x <= x_ref_r & z == z0)   = [];

%==========================================================================
% DEFINE PFIX REFINED REGION
%==========================================================================
pfix_vertex = [x_ref_l  z0 - d_ref ; ...
               x_ref_r  z0 - d_ref ; ...
               x_ref_r  z0         ; ...
               x_ref_l  z0 ]; % vertex points COUNTER-CLOCKWISE
nfix_vertex = size(pfix_vertex,1);
nseg        = nfix_vertex;

%==========================================================================
% DEFINE EACH BOUNDARY SEGMENT
%==========================================================================
seg_connect    = [(1:nseg)',[(2:nseg)';1]];       % This is the connectivity matrix for the fixed points
segBC          = zeros(nseg,6);                   % pre-allocate memory
segBC(:,[1 2]) = pfix_vertex(seg_connect(:,1),:); % p0 = (x0,y0) for seg_start
segBC(:,[3 4]) = pfix_vertex(seg_connect(:,2),:); % p1 = (x1,y1) for seg_end
segBC(:,[5 6]) = GetLineSlope(segBC(:,[1 2]),segBC(:,[3 4]),nseg); % slope_info = [a,b]

%==========================================================================
% COMPUTE PFIX (only valid for rectangular embedded high resolution region)
%==========================================================================
pfix = [];
for i = 1:nseg
    pfix_this_segment_ref = regular_bnd_nodes(segBC,l0_ref,i);
    pfix                  = [pfix; pfix_this_segment_ref];
end

%==========================================================================
% OUTPUT DATA
%==========================================================================
MESH.pfix_mesh_vertex = MESH.pfix;
MESH.pfix_ref_vertex  = pfix_vertex;
MESH.pfix             = [MESH.pfix; pfix_vertex; pfix];
MESH.pbnd             = pbnd;
MESH.pseg             = pseg;
MESH.GCOORD           = [MESH.pfix;MESH.pbnd;MESH.pint];
MESH.EL2NOD           = delaunay(MESH.GCOORD);
MESH.q                = mesh_quality(MESH.GCOORD,MESH.EL2NOD);

% SETTINGS.show_figs = 1;
% SETTINGS.save_figs = 0;
% plot_first_guess_mesh(MESH,SETTINGS)

end % END OF FUNCTION pfix_refined_region