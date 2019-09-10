function MESH = mesh_info(MESH,GUIDE_MESH)
% Usage: MESH = mesh_info(MESH,GUIDE_MESH)
%
% Purpose:
%   Compute number of nodes and elements inside the different regions
%   (coarse, transition and refined) for a 10 nodel mesh
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
x0         = GUIDE_MESH.x0;          % point around which the refined and transition zones are defined
z0         = GUIDE_MESH.z0;          % point around which the refined and transition zones are defined
d_tran     = GUIDE_MESH.d_tran;      % transition zone depth (km)
l_tran     = GUIDE_MESH.l_tran;      % length of transition zone (km)
x_tran_l   = x0 - l_tran/2;          % left boundary of the transition zone
x_tran_r   = x0 + l_tran/2;          % right boundary of the transition zone
d_ref      = GUIDE_MESH.d_ref;       % refined zone depth (km)
l_ref      = GUIDE_MESH.l_ref;       % length of refined zone (km)
x_ref_l    = x0 - l_ref/2;           % left boundary of the refined zone
x_ref_r    = x0 + l_ref/2;           % right boundary of the refined zone

%==========================================================================
% CREATE A 6 NODEL MESH FROM A 3 NODEL MESH
%==========================================================================
EL2NOD            = uint32(MESH.EL2NOD');
GCOORD            = MESH.GCOORD';
PointID           = point_id_rectbox(GCOORD);
[GCOORD,EL2NOD,~] = trimesh_p1_to_p2(GCOORD,EL2NOD,PointID);
GCOORD            = GCOORD';
EL2NOD            = EL2NOD';

%==========================================================================
% COMPUTE EACH ELEMENT BARYCENTER (IT DETERMINES THE ELEMENT POSITION)
%==========================================================================
x_bary            = (GCOORD(EL2NOD(:,1),1) + GCOORD(EL2NOD(:,2),1) + GCOORD(EL2NOD(:,3),1))/3;
z_bary            = (GCOORD(EL2NOD(:,1),2) + GCOORD(EL2NOD(:,2),2) + GCOORD(EL2NOD(:,3),2))/3;
GCOORD_bary       = [x_bary z_bary];

%==========================================================================
% COMPUTE STATISTICS FOR COARSE REGION
%==========================================================================
EL2NOD_coarse     = EL2NOD;
% Remove those elements that are inside the transition zone
EL2NOD_coarse(x_bary > x_tran_l & x_bary < x_tran_r & ...
              z_bary > z0 - d_tran                  ,:) = [];
% Select those nodes belonging to elements in the coarse region
GCOORD_coarse     = unique(GCOORD(EL2NOD_coarse,:),'rows','stable');
% Remove those nodes that are outside of the coarse region (even if the
% element they belong is considered inside the coarse region)
x_coarse          = GCOORD_coarse(:,1);
z_coarse          = GCOORD_coarse(:,2);
GCOORD_coarse(x_coarse > x_tran_l & x_coarse < x_tran_r & ...
              z_coarse > z0 - d_tran                    ,:) = [];

%==========================================================================
% COMPUTE STATISTICS FOR TRANSITION REGION
%==========================================================================
EL2NOD_tran       = EL2NOD;
% Take those elements that are inside the transition zone
EL2NOD_tran       = EL2NOD_tran(x_bary >= x_tran_l  & x_bary <= x_tran_r & ...
                                z_bary >= z0 - d_tran,:);
GCOORD_bary_tran  = GCOORD_bary;
GCOORD_bary_tran  = GCOORD_bary_tran(x_bary >= x_tran_l  & x_bary <= x_tran_r & ...
                                     z_bary >= z0 - d_tran,:);
% Remove those elements that are inside the refined zone
x_tran            = GCOORD_bary_tran(:,1);
z_tran            = GCOORD_bary_tran(:,2);
EL2NOD_tran(x_tran > x_ref_l & x_tran < x_ref_r & ...
            z_tran > z0 - d_ref                    ,:) = [];
% Select those nodes belonging to elements in the transition region
GCOORD_tran       = unique(GCOORD(EL2NOD_tran,:),'rows','stable');
% Remove those nodes that are outside of the transition region (even if the
% element they belong is considered inside the transition region)
% Take those nodes that are inside the transition zone
x_tran_temp       = GCOORD_tran(:,1);
z_tran_temp       = GCOORD_tran(:,2);
GCOORD_tran       = GCOORD_tran(x_tran_temp >= x_tran_l  & x_tran_temp <= x_tran_r & ...
                           z_tran_temp >= z0 - d_tran,:);
% Remove those elements that are inside the refined zone
x_tran_temp       = GCOORD_tran(:,1);
z_tran_temp       = GCOORD_tran(:,2);
GCOORD_tran(x_tran_temp > x_ref_l & x_tran_temp < x_ref_r & ...
            z_tran_temp > z0 - d_ref                      ,:) = []; 

%==========================================================================
% COMPUTE STATISTICS FOR REFINED REGION
%==========================================================================
EL2NOD_ref        = EL2NOD;
% Take those elements that are inside the refined zone
EL2NOD_ref        = EL2NOD_ref(x_bary >= x_ref_l  & x_bary <= x_ref_r & ...
                              z_bary >= z0 - d_ref,:);
% Select those nodes belonging to elements in the refined region
GCOORD_ref        = unique(GCOORD(EL2NOD_ref,:),'rows','stable');
% Remove those nodes that are outside of the refined region (even if the
% element they belong is considered inside the refined region)
x_ref             = GCOORD_ref(:,1);
z_ref             = GCOORD_ref(:,2);
GCOORD_ref        = GCOORD_ref(x_ref >= x_ref_l  & x_ref <= x_ref_r & ...
                               z_ref >= z0 - d_ref,:);

%==========================================================================
% DISPLAY INFORMATION
%==========================================================================
% Percentage of nodes and elements inside the coarse region
pct_nodes_coarse    = (size(GCOORD_coarse,1)/(size(GCOORD,1)))*100;
pct_elements_coarse = (size(EL2NOD_coarse,1)/(size(EL2NOD,1)))*100;

% Percentage of nodes and elements inside the transition region
pct_nodes_tran      = (size(GCOORD_tran,1)/(size(GCOORD,1)))*100;
pct_elements_tran   = (size(EL2NOD_tran,1)/(size(EL2NOD,1)))*100;

% Percentage of nodes and elements inside the refined region
pct_nodes_ref       = (size(GCOORD_ref,1)/(size(GCOORD,1)))*100;
pct_elements_ref    = (size(EL2NOD_ref,1)/(size(EL2NOD,1)))*100;

MESH_INFO.coarse = [pct_nodes_coarse pct_elements_coarse];
MESH_INFO.tran   = [pct_nodes_tran   pct_elements_tran  ];
MESH_INFO.ref    = [pct_nodes_ref    pct_elements_ref   ];

display_mesh_info(MESH_INFO)

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pct_nodes_coarse    = pct_nodes_coarse;
MESH.pct_elements_coarse = pct_elements_coarse;
MESH.pct_nodes_tran      = pct_nodes_tran;
MESH.pct_elements_tran   = pct_elements_tran;
MESH.pct_nodes_ref       = pct_nodes_ref;
MESH.pct_elements_ref    = pct_elements_ref;

end % END OF FUNCTION mesh_info

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function display_mesh_info(MESH_INFO)
% Usage: display_mesh_info(MESH_INFO)
% 
% Purpose: 
%   Displays number of nodes and elements inside the different regions
%   (coarse, transition and refined)  in structure "MESH_INFO" in Matlab's
%   terminal.
%
% Input:
%   MESH_INFO : [structure] : number of nodes and elements inside the
%                             different regions (coarse, transition and
%                             refined)
%
% Output:
%   none (Output only in Matlab terminal)
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JMT Jul 2016: Addapted to MESH_2D_SPRING

varnames = fieldnames(MESH_INFO);
nvar     = length(varnames);
fprintf('\n\n  INFO for a 6 nodel mesh\n');
fprintf(' --------------------------------\n');
fprintf(' Region  | %% nodes | %% elements |\n');
fprintf(' --------------------------------\n');
for i=1:nvar
    varname = varnames{i};
    l       = 8 - length(varname);
    m       = 10 - length(num2str(MESH_INFO.(varname)(1)));
    n       = 11 - length(num2str(MESH_INFO.(varname)(2)));
    fprintf(' %s%s| ',varname,repmat(' ',1,l));
    fprintf(' %4.1f%s|    %4.1f%s|\n',MESH_INFO.(varname)(1),repmat(' ',1,m),MESH_INFO.(varname)(2),repmat(' ',1,n));
end
fprintf(' --------------------------------');

% save the information
fid = fopen('INFO_6nodel_mesh.txt','w');
fprintf(fid,' --------------------------------\n');
fprintf(fid,' Region  | %% nodes | %% elements |\n');
fprintf(fid,' --------------------------------\n');
for i=1:nvar
    varname = varnames{i};
    l       = 8 - length(varname);
    m       = 10 - length(num2str(MESH_INFO.(varname)(1)));
    n       = 11 - length(num2str(MESH_INFO.(varname)(2)));
    fprintf(fid,' %s%s| ',varname,repmat(' ',1,l));
    fprintf(fid,' %4.1f%s|    %4.1f%s|\n',MESH_INFO.(varname)(1),repmat(' ',1,m),MESH_INFO.(varname)(2),repmat(' ',1,n));
end
fclose(fid);

end % END OF FUNCTION display_mesh_info