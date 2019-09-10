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
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% GENERAL SETTINGS
%==========================================================================
% MESH DIMENSIONS
SETTINGS.depth      = 2900;  % domain depth (km)
SETTINGS.length     = 40000; % domain length (km)
% Point around which the domain is created                
SETTINGS.x0         = 0;          % (km)
SETTINGS.z0         = 0;          % (km)

% TOLERANCES
SETTINGS.itmax      = 20;   % iteration tolerance to end the loop
SETTINGS.q_tol      = 0.45; % quality tolerance for the worst triangle
SETTINGS.q_smooth   = 0.5;    % quality tolerance for smoothing interior nodes
SETTINGS.q_balloon  = 0.5;  % quality tolerance for balloon forces
SETTINGS.mean_q_tol = 0.89;  % quality tolerance for q mean
SETTINGS.mean_misfit_bar_length_tol = 0.025; % tolerance for the misfit of springs (it sets
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
                GUIDE_MESH.x0        = 0; % (km)
                GUIDE_MESH.z0        = 0; % (km)
                
                % Coarse zone
                GUIDE_MESH.l0_coarse = 1500;   % desired spring length (km) for coarse zone
                
                % Transition zone
                GUIDE_MESH.d_tran    = 2900;   % transition zone depth (km)
                GUIDE_MESH.l_tran    = 8000;   % transition zone length (km)
                
                % Refined zone
                GUIDE_MESH.l0_ref    = 7.5;     % desired spring length (km) for refined zone
                GUIDE_MESH.d_ref     = 300;    % refined zone depth (km)
                GUIDE_MESH.l_ref     = 3333.3; % refined zone length (km)
            case 'load'
                % Full path to file containing the guide mesh is required!
                SETTINGS.guide_mesh_file = 'SETUP_TEST/GUIDE_MESH';
                GUIDE_MESH               = struct([]); % empty structure
            otherwise
                error('typo in SETTINGS.guide_mesh') 
        end
        
    case 'regular'
        % Desired spring length for a regular mesh
        SETTINGS.h0   = 1500;       % (km)  
        
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

if strcmp(SETTINGS.first_guess,'load')
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
SETTINGS.outdir = [path2data 'RECT_MESH/Trash_01']; % full path to output folder

end % END OF FUNCTION mesh_parameters