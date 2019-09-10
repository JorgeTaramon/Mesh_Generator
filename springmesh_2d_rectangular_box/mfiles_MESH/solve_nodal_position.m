function MESH = solve_nodal_position(MESH,GUIDE_MESH,SETTINGS)
% Usage: MESH = solve_nodal_position(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Compute new nodal position after applying the Finite Element Method to
%   solve Hooke's law for all the springs of the mesh.
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
% ASSEMBLY STIFFNESS MATRIX AND FORCE VECTOR
%==========================================================================
[MESH,STIFFNESS,F_le] = assembly_stiffness(MESH,GUIDE_MESH,SETTINGS);

%==========================================================================
% CROSS BARS AND BALLOON FORCES
%==========================================================================
STIFFNESS.KK1_CrossBars = [];
STIFFNESS.KKi_CrossBars = [];
STIFFNESS.KKj_CrossBars = [];

if SETTINGS.cross_bars
    [STIFFNESS.KK1_CrossBars,STIFFNESS.KKi_CrossBars,STIFFNESS.KKj_CrossBars,F_le] = cross_bars(MESH,F_le);
end

if SETTINGS.balloon_forces
    [q_sorted,I]  = sort(MESH.q); % q sorted starting from the worst q
    EL2NOD_sorted = MESH.EL2NOD(I(q_sorted < SETTINGS.q_balloon),:); % select those elements which q < q_balloon
    if ~isempty(EL2NOD_sorted)
        F_balloon = balloon_forces(MESH,EL2NOD_sorted,GUIDE_MESH,SETTINGS);
        F_le      = F_le - F_balloon; % we have to subtract F_balloon since F_le is in the rigth  
                                        % hand side and F_balloon is originaly in the left hand side
    end
end

%==========================================================================
% BC METHOD AND SOLVE FOR UNKNOWNS (NODAL POSITIONS)
%==========================================================================
switch SETTINGS.bc_method
    case 'projection'
        [MESH] = projection(MESH,STIFFNESS,F_le,SETTINGS);
    case 'uzawa'
        [MESH] = uzawa(MESH,STIFFNESS,F_le,SETTINGS); % need to be coded!!!
    otherwise
        error('typo in SETTINGS.bc_method, see mesh_parameters.m')
end

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
[MESH.L0,MESH.L,MESH.bars,~,~] = bar_length(MESH,GUIDE_MESH,SETTINGS);
MESH.rel_change                = (MESH.L - MESH.L0)./MESH.L0; % relative bar-length change
MESH.rel_change_abs            = abs(MESH.rel_change);        % absolute value of relative bar-length change
MESH.mean_misfit_bar_length    = sum(MESH.rel_change_abs)/size(MESH.L,1);
MESH.q                         = mesh_quality(MESH.GCOORD,MESH.EL2NOD);
MESH.rho0                      = sqrt(2)./(MESH.L0.^2); % desired nodal desnsity
MESH.rho                       = sqrt(2)./(MESH.L.^2);  % actual nodal desnsity

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.save_figs || SETTINGS.show_figs
    % Plot a histogram of percentage of misfit for the bars
    figure(4); clf
    plot_misfit_bar_length_to_check_nodal_density(MESH,SETTINGS)
end

end %END OF FUNCTION solve_nodal_position