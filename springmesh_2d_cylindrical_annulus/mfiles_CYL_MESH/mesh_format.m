function mesh_format(MESH,GUIDE_MESH,SETTINGS)
% Usage: mesh_format(MESH,GUIDE_MESH,SETTINGS)
% 
% Purpose: 
%   Make a mesh created by MESH_2D_SPRING_CYL (3 nodel) compatible with the
%   2D Finite Element convection code in polar coordinates M2TRI_CYL and
%   save it (6 nodel)
%
% Input:
%   MESH : [structure] : FE mesh in MESH_2D_SPRING_CYL format
%
% Output:
%   none 
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% CREATE VARIABLES FROM A 3 NODEL MESH
%==========================================================================
EL2NOD          = uint32(MESH.EL2NOD');
GCOORD          = MESH.GCOORD';
GCOORD_POL      = cartesian2polar_rad(GCOORD);
r_int           = SETTINGS.r_int;
r_ext           = SETTINGS.r_ext;
PointID         = zeros(1,size(GCOORD,2));
PointID(sqrt(GCOORD(1,:).^2+GCOORD(2,:).^2) < r_int+0.0001) = 201;
PointID(sqrt(GCOORD(1,:).^2+GCOORD(2,:).^2) > r_ext-0.0001) = 203;
% PhaseID      = 1;

%==========================================================================
% CREATE A 6 NODEL MESH WITH CURVED EDGES FROM 3 NODEL MESH
%==========================================================================
[GCOORD,~,EL2NOD,PointID] = ...
    trimesh_p1_to_p2_cyl(GCOORD,GCOORD_POL,EL2NOD,PointID);

nnod     = size(GCOORD,2);
name     = strcat('EarthCyl_n',num2str(nnod));
filename = [SETTINGS.outdir '/' name];

%==========================================================================
% CREATE THE STRUCTURE FOR A 6 NODEL MESH
%==========================================================================
clear MESH
MESH.EL2NOD  = EL2NOD;
MESH.GCOORD  = GCOORD;
MESH.PointID = PointID;
PhaseID      = ones(1,size(EL2NOD,2));
% x_bary       = (GCOORD(1,EL2NOD(1,:)) + GCOORD(1,EL2NOD(2,:)) + GCOORD(1,EL2NOD(3,:)))/3;
% z_bary       = (GCOORD(2,EL2NOD(1,:)) + GCOORD(2,EL2NOD(2,:)) + GCOORD(2,EL2NOD(3,:)))/3;
% PhaseID(sqrt(x_bary.^2 + z_bary.^2) <= (4921)) = 2;
MESH.PhaseID = PhaseID;
MESH.name    = name;
if strcmp(SETTINGS.refinement,'guide_mesh')
    % SAVE DIMENSIONS OF REFINED ZONE
    MESH.d_ref  = GUIDE_MESH.d_ref; % depth
    MESH.w_ref  = GUIDE_MESH.w_ref; % width
    % SAVE ELEMENT LENGTH INSIDE REFINED ZONE
    MESH.l0_ref = GUIDE_MESH.l0_ref;
end
MESH.r_cmb  = r_int;
MESH.r_surf = r_ext;
save(filename,'MESH');

end % END OF FUNCTION mesh_format