function create_multigrid_surface_nodes()
% Usage: create_multigrid_surface_nodes()
% 
% Purpose: 
%   Recursively refine the mesh by splitting each tetrahedron into 8 elements
%   
% Input:
%
% Output:
%   Surface nodes saved in txt file
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

pdir = pwd;
cd('..');
addpath([pwd '/mfiles_SPH_MESH']);
cd(pdir);

%==========================================================================
% LOAD MESH DATA
%==========================================================================

% Specify output directory for the mesh(folder where data is located)
outdir_mesh   = ['/Users/jorge/Tests/SPH_MESH/Trash_00' filesep];
filename_mesh = 'EarthShell_n37092.mat';

fprintf(' Loading mesh data and settings...');
try
    data  = load([outdir_mesh filename_mesh]);
catch
    error(' Could no open file "%s"\n',[outdir_mesh 'MESH.mat']);
end
MESH = data.MESH;
clear data

try
    data  = load([outdir_mesh 'GUIDE_MESH.mat']);
catch
    error(' Could no open file "%s"\n',[outdir_mesh 'GUIDE_MESH.mat']);
end
GUIDE_MESH = data.GUIDE_MESH;
clear data

try
    data  = load([outdir_mesh 'SETTINGS.mat']);
catch
    error(' Could no open file "%s"\n',[outdir_mesh 'SETTINGS.mat']);
end
SETTINGS = data.SETTINGS;
clear data2
fprintf(' done\n');

GCOORD_c            = MESH.GCOORD;
GCOORD_SPH_c        = cartesian2spherical_rad(GCOORD_c);
EL2NOD_c            = MESH.EL2NOD;
PointID_c           = MESH.PointID;
PhaseID_c           = ones(1,size(EL2NOD_c,2),'int32');
DB_indices{1}       = 301;
DB_indices{2}       = 306;
r_int               = MESH.r_cmb;
r_ext               = MESH.r_surf;
SETTINGS.theta_cone = 45;
SETTINGS.show_figs  = 0;
nnod                = size(GCOORD_c,2);
prompt              = 'Type number of multigrid levels (integer): ';
nmg                 = input(prompt);

%==========================================================================
% REFINE THE MESH BY SPLITTING EACH TETRAHEDRON RECURSIVELY
%==========================================================================
for img=nmg:-1:2
    
    % CALCULATE REFINED MESH
    [GCOORD,GCOORD_SPH,EL2NOD,PointID,PhaseID] = ...
        tetmesh_refine_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices,r_int,r_ext,SETTINGS);

    % SAVE CURRENT MESH AS 'COARSE'
    if img>2
        GCOORD_c     = GCOORD;
        GCOORD_SPH_c = GCOORD_SPH;
        EL2NOD_c     = EL2NOD;
        PointID_c    = PointID;
        PhaseID_c    = PhaseID;
    end
end

%==========================================================================
% CREATE THE STRUCTURE FOR A 10 NODEL MESH WITH CURVED EDGES
%==========================================================================
clear MESH
MESH.EL2NOD  = EL2NOD;
MESH.GCOORD  = GCOORD;
MESH.PointID = PointID;
% MESH.PhaseID = PhaseID;
if strcmp(SETTINGS.refinement,'guide_mesh')
    % SAVE DIMENSIONS OF REFINED ZONE
    MESH.d_ref        = GUIDE_MESH.d_ref; % depth
    MESH.w_ref        = GUIDE_MESH.w_ref; % width (North-South)
    MESH.l_ref        = GUIDE_MESH.l_ref; % length (East-West)
    % SAVE ELEMENT LENGTH INSIDE REFINED ZONE
    MESH.l0_ref       = GUIDE_MESH.l0_ref;
    % SAVE DIMENSIONS OF TRANSITION ZONE
    MESH.d_tran       = GUIDE_MESH.d_tran; % depth
    MESH.w_tran       = GUIDE_MESH.w_tran; % width (North-South)
    MESH.l_tran       = GUIDE_MESH.l_tran; % length (East-West)
    % SAVE POINT AROUND WHICH REFINED AND TRANSITION REGIONS ARE CREATED
    MESH.theta0       = GUIDE_MESH.theta0; % colatitude (degrees) of the point around which the refined and transition zones are defined
    MESH.phi0         = GUIDE_MESH.phi0;   % longitude (degrees) of the point around which the refined and transition zones are defined
    % SAVE POINT AROUND WHICH THE VELOCITY BCs ARE TAKEN (e.g., South Atlantic MOR 130 Myr ago) AND FINITE ROTATION MATRICES 
    MESH.theta_center = GUIDE_MESH.theta_center;
    MESH.phi_center   = GUIDE_MESH.phi_center;
    MESH.RR1_center   = rotation_matrix_to_center_high_res_region(GUIDE_MESH.EP1_lat,GUIDE_MESH.EP1_lon,GUIDE_MESH.EP1_angle);
    MESH.RR2_center   = rotation_matrix_to_center_high_res_region(GUIDE_MESH.EP2_lat,GUIDE_MESH.EP2_lon,GUIDE_MESH.EP2_angle);
end
MESH.r_cmb  = r_int;
MESH.r_surf = r_ext;
save([outdir_mesh '/' strcat('EarthShell_n',num2str(nnod),'_nmg',num2str(nmg))],'MESH');

%==========================================================================
% CREATE .txt FILE FOR SURFACE NODES (latitude, longitude)
%==========================================================================
if strcmp(SETTINGS.refinement,'guide_mesh')
    GCOORD_centered = MESH.RR2_center * MESH.RR1_center * MESH.GCOORD;  % rotate points to be centered around the point of interest (theta_center,phi_center)
    GCOORD_SPH      = transpose(cartesian2spherical(GCOORD_centered')); % spherical coordinates (theta,phi,r)
else
    GCOORD_SPH = transpose(cartesian2spherical(MESH.GCOORD')); % spherical coordinates (theta,phi,r)
end
inod_surf      = find(MESH.PointID == 306);                    % nodes on the surface
lambda         = 90 - GCOORD_SPH(1,inod_surf);                 % lambda = 90° - theta (latitude = 90° - colatitude)
phi            = GCOORD_SPH(2,inod_surf);                      % phi (longitude)
latlon_surf    = [lambda; phi];
cd('..');
fid = fopen('surf_nodes.txt','w');
fprintf(fid,'%d %d \n',latlon_surf);
fclose(fid);

txt_file_source_path = [pwd '/surf_nodes.txt']; 
cd GPlates
movefile(txt_file_source_path)

end % END OF FUNCTION create_multigrid_meshes

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [GCOORD,GCOORD_SPH,EL2NOD,PointID,PhaseID] = ...
    tetmesh_refine_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,PhaseID_c,DB_indices,r_int,r_ext,SETTINGS) 
     
    % Check if the mesh has linear (4-node) or quadratic (10-node) tetrahedral
    % elements. If quadratic, split each element into 8 linear sub-elements. 
    % =========================================================================
    nnodel = size(EL2NOD_c,1);
    if nnodel==10
        [EL2NOD_c,PhaseID_c] = tetmesh_p2_to_p1(GCOORD_c,EL2NOD_c,PhaseID_c);
    end
    
    % =======================================================================================
    % CREATE A 10 NODE MESH WITH CURVED EDGE ELEMENTS
    % =======================================================================================
    [GCOORD,GCOORD_SPH,EL2NOD,PointID,~,~,~,~] = ...
        tetmesh_p1_to_p2_sph(GCOORD_c,GCOORD_SPH_c,EL2NOD_c,PointID_c,DB_indices,r_int,r_ext,SETTINGS);
    
    if nnodel==4
        % Re-connect elements and nodes by creating a linear (4-node) element
        % connectivity matrix (only if input was a linear connectivity)
        % =========================================================================
        [EL2NOD,PhaseID] = tetmesh_p2_to_p1(GCOORD_c,EL2NOD_c,PhaseID_c);
    else
        PhaseID          = PhaseID_c;
    end
end