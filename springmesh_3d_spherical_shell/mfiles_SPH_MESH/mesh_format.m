function mesh_format(MESH,GUIDE_MESH,SETTINGS)
% Usage: mesh_format(MESH,GUIDE_MESH,SETTINGS)
% 
% Purpose: 
%   Make the mesh created by MESH_3D_SPRING_SPH (4 nodel) compatible with 
%   the 3D Finite Element convection code in spherical coordinates 
%   M3TET_SPH and save it (10 nodel).
%   
% Input:
%   MESH       : [structure] : FE mesh in MESH_3D_SPRING_SPH format
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   none 
%
% JMT Jul 2016
% JMT Aug 2016: Now it saves surfaces nodes in a .txt file to get the plate
%               velocities from GPlates (before we need to transform the
%               .txt file into a .gpml file)
% JMT Oct 2017: Now the output quadratic mesh (10 nodel) can have straight
%               edges or curved edges
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% CREATE VARIABLES FROM A 4 NODEL MESH
%==========================================================================
EL2NOD            = uint32(MESH.EL2NOD');
GCOORD            = MESH.GCOORD';
GCOORD_SPH        = cartesian2spherical_rad(GCOORD);
r_int             = SETTINGS.r_int;
r_ext             = SETTINGS.r_ext;
% Check if node connectivity is ok; If a negative volume is calculated
% (when assuming Hughes' notation), convert the numbering to Hughes' notation 
el_vol            = calc_tetra_volume(GCOORD,EL2NOD);
iel               = el_vol<1;
if any(iel)
    EL2NOD(:,iel) = EL2NOD([1 2 4 3]',iel); % change connectivity
end
PointID           = zeros(1,size(GCOORD,2));
PointID(sqrt(sum(GCOORD.^2,1)) < r_int+0.0001) = 301; % point IDs for CMB
PointID(sqrt(sum(GCOORD.^2,1)) > r_ext-0.0001) = 306; % point IDs for Earth's surface
% PhaseID           = 1;
DB_indices{1}     = 301;
DB_indices{2}     = 306;

%==========================================================================
% CREATE A 10 NODEL MESH FROM 4 NODEL MESH
%==========================================================================
switch SETTINGS.edges_output_mesh
    case 'straight' % 10 nodel mesh with straight edges
        [GCOORD,EL2NOD,PointID] = tetmesh_p1_to_p2(GCOORD,EL2NOD,PointID,DB_indices);
    case 'curved'   % 10 nodel mesh with curved edges
        SETTINGS.theta_cone = 45;
        SETTINGS.show_figs  = 0;
        [GCOORD,~,EL2NOD,PointID,~,~,~,~] = ...
            tetmesh_p1_to_p2_sph(GCOORD,GCOORD_SPH,EL2NOD,PointID,DB_indices,r_int,r_ext,SETTINGS);
end
nnod                    = size(GCOORD,2);
name                    = strcat('EarthShell_n',num2str(nnod));
filename                = [SETTINGS.outdir '/' name];

%==========================================================================
% CREATE THE STRUCTURE FOR A 10 NODEL MESH WITH CURVED EDGES
%==========================================================================
clear MESH
MESH.EL2NOD  = EL2NOD;
MESH.GCOORD  = GCOORD;
MESH.PointID = PointID;
% MESH.PhaseID = PhaseID;
MESH.name    = name;
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
save(filename,'MESH');

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
fid = fopen('surf_nodes.txt','w');
fprintf(fid,'%d %d \n',latlon_surf);
fclose(fid);

% %==========================================================================
% % TRANSFORM .txt FILE INTO .gpml FILE AND SAVE IT IN outdir
% %==========================================================================
% txt_file_source_path = [pwd '/surf_nodes.txt']; 
% cd GPlates
% movefile(txt_file_source_path)
% setenv('DYLD_LIBRARY_PATH','');
% system('./txt2gpml.sh')
% % % IMPORTANT!!
% % % If you get the following error:
% % % dyld: warning, LC_RPATH @executable_path/../lib in //Users/jorge/Dropbox/GEOMAR/springmesh_3d/GPlates/pygplates/pygplates.so being ignored in restricted program because of @executable_path
% % % dyld: warning, LC_RPATH @executable_path/../../../../lib in //Users/jorge/Dropbox/GEOMAR/springmesh_3d/GPlates/pygplates/pygplates.so being ignored in restricted program because of @executable_path
% % % Traceback (most recent call last):
% % %   File "PyGPlates_MeshNodePointsFromLatLonFile.py", line 23, in <module>
% % %      import pygplates
% % % ImportError: dlopen(//Users/jorge/Dropbox/GEOMAR/springmesh_3d/GPlates/pygplates/pygplates.so, 2): Library not loaded: @rpath/libpython2.7.dylib
% % %   Referenced from: //Users/jorge/Dropbox/GEOMAR/springmesh_3d/GPlates/pygplates/pygplates.so
% % %   Reason: image not found
% % % It seems to be a compatibility issue between Matlab Engine for python 2.7.12 and MacOSX 10.12 (Sierra) 
% % % see: https://uk.mathworks.com/matlabcentral/answers/283580-corrupted-version-of-matlab-engine-for-python-on-macosx-10-11
% % % However, if we open the terminal and type ./txt2gpml.sh in the same folder it works perfectly
% 
% % gpml_file_source_path      = [pwd '/points.gpml'];
% % gpml_file_destination_path = [SETTINGS.outdir '/' name '.gpml'];
% % movefile(gpml_file_source_path,gpml_file_destination_path)
% cd ..

end % END OF FUNCTION mesh_format