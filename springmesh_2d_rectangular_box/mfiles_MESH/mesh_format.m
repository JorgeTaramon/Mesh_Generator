function mesh_format(MESH,GUIDE_MESH,SETTINGS)
% Usage: mesh_format(MESH,GUIDE_MESH,SETTINGS)
% 
% Purpose: 
%   Make a mesh created by MESH_2D_SPRING (3 nodel) compatible with the
%   2D Finite Element convection code in Cartesian coordinates M2TRI and
%   save it (6 nodel)
%
% Input:
%   MESH : [structure] : FE mesh in MESH_2D_SPRING_CYL format
%
% Output:
%   none 
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% CREATE VARIABLES FROM A 3 NODEL MESH
%==========================================================================
EL2NOD          = uint32(MESH.EL2NOD');
GCOORD          = MESH.GCOORD';
PointID         = point_id_rectbox(GCOORD);
% PhaseID      = 1;

%==========================================================================
% CREATE A 6 NODEL MESH FROM 3 NODEL MESH
%==========================================================================
[GCOORD,EL2NOD,~] = trimesh_p1_to_p2(GCOORD,EL2NOD,PointID);
nnod              = size(GCOORD,2);
name              = strcat('Rectbox_n',num2str(nnod));
filename          = [SETTINGS.outdir '/' name];

%==========================================================================
% CREATE THE STRUCTURE FOR A 6 NODEL MESH
%==========================================================================
clear MESH
MESH.EL2NOD  = EL2NOD;
MESH.GCOORD  = GCOORD;
MESH.PointID = PointID;
% MESH.PhaseID = PhaseID;
MESH.name    = name;
if strcmp(SETTINGS.refinement,'guide_mesh')
    % SAVE DIMENSIONS OF REFINED ZONE
    MESH.d_ref  = GUIDE_MESH.d_ref; % depth
    MESH.l_ref  = GUIDE_MESH.l_ref; % length
    % SAVE ELEMENT LENGTH INSIDE REFINED ZONE
    MESH.l0_ref = GUIDE_MESH.l0_ref;
end
save(filename,'MESH');

end % END OF FUNCTION mesh_format