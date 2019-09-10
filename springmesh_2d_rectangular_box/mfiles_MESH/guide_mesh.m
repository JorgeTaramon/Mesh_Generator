function GUIDE_MESH = guide_mesh(SETTINGS,GUIDE_MESH)
% Usage: GUIDE_MESH = guide_mesh(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Generation of a guide-mesh for a 2D rectangular mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   GUIDE_MESH : [structure] : structure containing guide mesh settings and
%                              the guide mesh itself
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth            = SETTINGS.depth;         % domain depth (km)
length           = SETTINGS.length;        % domain length (km)
x0               = GUIDE_MESH.x0;          % point around which the refined and transition zones are defined
z0               = GUIDE_MESH.z0;          % point around which the refined and transition zones are defined
l0_coarse        = GUIDE_MESH.l0_coarse;   % desired length (km) for coarse zone (constant)
l0_tran          = GUIDE_MESH.l0_coarse;   % desired length (km) for the boundary of transition zone
l0_ref           = GUIDE_MESH.l0_ref;      % desired length (km) for refined zone (constant)
% Transition Zone
l_tran           = GUIDE_MESH.l_tran;      % length of transition zone (km)
x_tran_l         = x0 - l_tran/2;          % left boundary of the transition zone
x_tran_r         = x0 + l_tran/2;          % right boundary of the transition zone
% Refined Zone
d_ref            = GUIDE_MESH.d_ref;       % max depth in the refined zone (km)
l_ref            = GUIDE_MESH.l_ref;       % length of refined zone (km)
x_ref_l          = x0 - l_ref/2;           % left boundary of the refined zone
x_ref_r          = x0 + l_ref/2;           % right boundary of the refined zone

%==========================================================================
% GENERATE THE GUIDE-MESH
%==========================================================================
% Define points at a certain level (CMB) using the coordinates for transition and refined zones
x1               = [x0 - length/2; x_tran_l; x_ref_l; x_ref_r; x_tran_r; x0 + length/2];
z1               = repmat(z0 - depth,size(x1));
p_CMB_level      = [x1 z1];

% Coarse zone
p_ref_level      = [p_CMB_level(:,1) (z0 - d_ref)*ones(size(p_CMB_level,1),1)];
p_surf_level     = [p_CMB_level(:,1)           z0*ones(size(p_CMB_level,1),1)];
p_coarse         = [p_CMB_level; p_ref_level; p_surf_level];
p_coarse(p_coarse(:,1) >= x_tran_l &...
         p_coarse(:,1) <= x_tran_r,:) = []; % remove those nodes in transition and refined zones

% Transition zone 
p_CMB_level      = p_CMB_level(p_CMB_level(:,1) >= x_tran_l & p_CMB_level(:,1) <= x_tran_r,:);
p_ref_level      = [p_CMB_level(:,1) (z0 - d_ref)*ones(size(p_CMB_level,1),1)];
p_surf_level     = [p_CMB_level(:,1)           z0*ones(size(p_CMB_level,1),1)];
p_tran           = [p_CMB_level; p_ref_level; p_surf_level];
p_tran(p_tran(:,1) >= x_ref_l & p_tran(:,1) <= x_ref_r & ...
       p_tran(:,2) >= z0 - d_ref,:) = []; % remove those nodes in refined zone

% Refined zone 
p_CMB_level      = p_CMB_level(p_CMB_level(:,1) >= x_ref_l & p_CMB_level(:,1) <= x_ref_r,:);
p_ref_level      = [p_CMB_level(:,1) (z0 - d_ref)*ones(size(p_CMB_level,1),1)];
p_surf_level     = [p_CMB_level(:,1)           z0*ones(size(p_CMB_level,1),1)];
p_ref       	 = [p_ref_level; p_surf_level];

% Create GCOORD_GUIDE and EL2NOD_GUIDE
GCOORD_GUIDE     = [p_coarse;p_tran;p_ref];              % guide mesh nodes in polar coordinates (theta,r)
[GCOORD_GUIDE,I] = unique(GCOORD_GUIDE,'rows','stable'); % remove repeated nodes
EL2NOD_GUIDE     = delaunay(GCOORD_GUIDE);               % create connectivity matrix

%==========================================================================
% SET THE DESIRED LENGTH (L0) FOR EVERY NODE OF THE GUIDE-MESH
%==========================================================================
L0_guide_coarse  = l0_coarse*ones(size(p_coarse,1),1); % desired bar length for the position of each node in the coarse zone
L0_guide_tran    = l0_tran*ones(size(p_tran,1),1);     % desired bar length for the position of each node in the transition zone
L0_guide_ref     = l0_ref*ones(size(p_ref,1),1);       % desired bar length for the position of each node in the refined zone
L0_guide         = [L0_guide_coarse; L0_guide_tran; L0_guide_ref];
L0_guide         = L0_guide(I); % take only the L0 values for the right number of nodes(after removing repeated nodes)

%==========================================================================
% CREATE OUTPUT STRUCTURE
%==========================================================================
GUIDE_MESH.GCOORD_GUIDE   = GCOORD_GUIDE;
GUIDE_MESH.EL2NOD_GUIDE   = EL2NOD_GUIDE;
GUIDE_MESH.L0_guide       = L0_guide;
GUIDE_MESH.p_coarse_guide = p_coarse;
GUIDE_MESH.p_tran_guide   = p_tran;
GUIDE_MESH.p_ref_guide    = p_ref;
GUIDE_MESH.x0             = x0;
GUIDE_MESH.z0             = z0;
GUIDE_MESH.depth          = depth;
GUIDE_MESH.length         = length;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    plot_guide_mesh(GUIDE_MESH,SETTINGS)
end
end % END OF FUNCTION guide_mesh