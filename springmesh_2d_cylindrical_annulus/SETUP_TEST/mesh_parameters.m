function [SETTINGS,GUIDE_MESH] = mesh_parameters()
% Usage: [SETTINGS,GUIDE_MESH] = mesh_parameters()
%
% Purpose:
%   Define mesh and guide mesh parameters
%
% Input:
%   none
%
% Output:
%   SETTINGS    : [structure] : structure containing mesh settings
%   GUIDE_MESH  : [structure] : structure containing guide mesh settings
%
% JMT May 2016
% JMT Jun 2017: cleaned up
% JMT Oct 2017: added option to create half a cylindrical annulus mesh for
%               axisymmetric code
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% GENERAL SETTINGS
%==========================================================================
% MESH
SETTINGS.r_ext      = 6371; % outer radius (km) of cylindrical annulus mesh
SETTINGS.r_int      = 3471; % inner radius (km) of cylindrical annulus mesh. 
                            % Set it to 0 to create a circular mesh without
                            % the inner hole (only for regular mesh)
SETTINGS.mesh       = 'full'; % set the output mesh
                              % 'full'   -> mesh compatible with M2TRI_CYL
                              % 'axisym' -> mesh compatible with axisymmetric version of M2TRI_CYL
% TOLERANCES
SETTINGS.itmax      = 20;   % iteration tolerance to end the loop
SETTINGS.q_tol      = 0.40; % quality tolerance for the worst triangle
SETTINGS.q_smooth   = 0.5;  % quality tolerance for smoothing interior nodes
SETTINGS.q_balloon  = 0.5;  % quality tolerance for balloon forces
SETTINGS.mean_q_tol = 0.93; % quality tolerance for q mean
SETTINGS.mean_misfit_bar_length_tol = 0.04; % tolerance for the misfit of springs (it sets
                                            % when the right nodal density is achieved) 
% SOLVERS
SETTINGS.cross_bars     = 0; % whether to use cross bars (1--yes, 0--no)
                             % steadier structure in each triangle.
SETTINGS.balloon_forces = 0; % whether to use ballon-forces (1--yes, 0--no)
                             % (for having better triangles)
SETTINGS.bc_method      = 'projection'; % keep projection!! uzawa not coded 
                                        % 'projection' --> fixing bnd nodes on tangent line 
                                        %                  to circumference at each node 
                                        % 'uzawa'      --> penalty method   

%==========================================================================
% REFINEMENT SETTINGS
%==========================================================================
SETTINGS.refinement = 'guide_mesh'; % set refinement type
                                    % 'guide_mesh' --> refine the mesh using a guide-mesh
                                    % 'regular'    --> create a regular mesh
switch SETTINGS.refinement
    case 'guide_mesh'
        SETTINGS.guide_mesh = 'new'; % set either create or load a guide-mesh
                                     % 'new'  --> create a new guide mesh
                                     % 'load' --> load a pre-existing guide mesh
        SETTINGS.bnds_at_ref_region = 'fixed'; % set whether creating boundaries for the refined region.
                                               % 'fixed' -> Nodes confined in those  boundaries are fixed
                                               % 'no'    -> No boundaries for refined region
        switch SETTINGS.guide_mesh
            case 'new'
                % Point around which refined and transition zones are created
                GUIDE_MESH.theta0    = 90;     % colatitude (degrees)
                GUIDE_MESH.r0        = 6371;   % radial distance (km)
                
                % Coarse zone
                GUIDE_MESH.l0_coarse = 2000;   % desired spring length (km) for coarse zone
                
                % Transition zone
                GUIDE_MESH.d_tran    = 2900;   % transition zone depth (km)
                GUIDE_MESH.w_tran    = 8000;   % transition zone width (km)
                
                % Refined zone
                GUIDE_MESH.l0_ref    = 10;     % desired spring length (km) for refined zone
                GUIDE_MESH.d_ref     = 300;    % refined zone depth (km)
                GUIDE_MESH.w_ref     = 3333.3; % refined zone width (km)
            case 'load'
                % Full path to file containing the guide mesh is required!
                SETTINGS.guide_mesh_file = 'SETUP_TEST/GUIDE_MESH';
                GUIDE_MESH               = struct([]); % empty structure
            otherwise
                error('typo in SETTINGS.guide_mesh') 
        end
        
    case 'regular'
        SETTINGS.h0   = 2000;       % desired spring length (km) for a regular mesh
        GUIDE_MESH    = struct([]); % empty structure
        
    otherwise
        error('typo in SETTINGS.refinement')
end

%==========================================================================
% FIRST GUESS SETTINGS
%==========================================================================
SETTINGS.first_guess = 'load';
    % 'new'          --> create a new first guess for the mesh
    % 'load'         --> load a pre-existing mesh
    % 'load_as_pfix' --> load a pre-existing mesh as pfix and create the rest 
    %                    of nodes according to the parameters chosen in SETTINGS
    %                    (useful to create a spherical shell with different layers)

if strcmp(SETTINGS.first_guess,'load') || ...
   strcmp(SETTINGS.first_guess,'load_as_pfix')
    % Load a first guess mesh from file.
    % Full path to file containing the guide mesh is required!
    % File must contain the following fields:
    % "pfix"  : coord of fixed nodes
    % "pbnd1" : coord of boundary nodes for inner boundary
    % "pbnd2" : coord of boundary nodes for outer boundary
    % "pint"  : coord of interior nodes
    SETTINGS.first_guess_file = 'SETUP_TEST/MESH';
end

%==========================================================================
% FIGURE SETTINGS
%==========================================================================
SETTINGS.iplot       = 1;     % plot file number to enter the main loop
SETTINGS.save_figs   = 0;     % 0--> figures not saved (faster)
                              % 1--> figures saved in output directory
SETTINGS.fig_type    = 'png'; % 'png' or 'fig'; format of saved figure
SETTINGS.show_figs   = 0;     % 0--> figures not shown on screen (faster)
                              % 1--> figures shown on screen

%==========================================================================
% OUTPUT FOLDER
%==========================================================================
path2data       = data_storage_2d(); % edit file "data_storage_2d" to add 
                                     % location on new computer
SETTINGS.outdir = [path2data 'CYL_MESH/Trash_01']; % full path to output folder

end % END OF FUNCTION mesh_parameters