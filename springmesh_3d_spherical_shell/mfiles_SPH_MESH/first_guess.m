function MESH = first_guess(SETTINGS,GUIDE_MESH)
% Usage: MESH = first_guess(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Define the first guess for the nodes of a 3D spherical shell mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   MESH       : [structure] : structure containing 1st guess for the mesh
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% CREATE BOUNDARY AND INTERIOR NODES
%==========================================================================
if ~isempty(GUIDE_MESH) == 1
    [pfix,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH);
    pint               = interior_nodes(SETTINGS,GUIDE_MESH);
else
    [pfix,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS);
    pint               = interior_nodes_no_guide(SETTINGS);
end

%==========================================================================
% OUTPUT DATA
%==========================================================================
GCOORD     = [pfix;pbnd1;pbnd2;pint];
EL2NOD     = delaunay(GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH = cartesian2spherical(GCOORD);
EL2NOD     = EL2NOD(~(sum(ismember(EL2NOD,find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8)),2)==4),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_SPH(:,3)-SETTINGS.r_int) < 1e-8);
% elements_with_4_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 4;
% remove_elements_with_4_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_4_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------
MESH.pfix   = pfix;
MESH.pbnd1  = pbnd1;
MESH.pbnd2  = pbnd2;
MESH.pint   = pint;
MESH.GCOORD = GCOORD;
MESH.EL2NOD = EL2NOD;
MESH.q      = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
MESH.s      = shape_measure(MESH.GCOORD,MESH.EL2NOD);
MESH.r_int  = SETTINGS.r_int;
MESH.r_ext  = SETTINGS.r_ext;

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    plot_first_guess_mesh(MESH,SETTINGS)
    MESH.iter      = 0;
    SETTINGS.iplot = 0;
    figure(15)
    plot_mesh(MESH,SETTINGS)
    figure(16)
    MESH.worst_q = min(MESH.q);
    MESH.mean_q  = mean(MESH.q);
    plot_statistic_log_v2(MESH,SETTINGS)
end

end % END OF FUNCTION first_guess