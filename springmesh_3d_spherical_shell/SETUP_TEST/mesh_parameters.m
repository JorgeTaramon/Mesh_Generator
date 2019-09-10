function [SETTINGS,GUIDE_MESH,INTERFACE] = mesh_parameters()
% Usage: [SETTINGS,GUIDE_MESH,INTERFACE] = mesh_parameters()
%
% Purpose:
%   Define mesh and guide mesh parameters
%
% Input:
%   none
%
% Output:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%   INTERFACE  : [structure] : structure containing interface settings
%
% JMT Jul 2015
% JMT Jun 2016: cleaned up
% JMT Jan 2017: reorganized and added the option to refine interface layers
% JMT Aug 2017: Added option to create internal boundaries for the refined
%               region. The nodes of the internal boundaries for the 
%               refined region can be 'fixed' or 'free' to move along the
%               internal boundaries.
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% GENERAL SETTINGS
%==========================================================================
% MESH DIMENSIONS
SETTINGS.r_ext      = 6371; % outer radius (km) of spherical shell
SETTINGS.r_int      = 3471; % inner radius (km) of spherical shell. 
                            % Set it to 0 to create a sphere instead of a shell
                            % (only for regular mesh)

% TOLERANCES
SETTINGS.itmax      = 20;   % iteration tolerance to end the loop
SETTINGS.q_tol      = 0.23; % quality tolerance for the worst tetrahedron
SETTINGS.q_bad      = 0.23; % quality tolerance for bad tetrahedrons
SETTINGS.q_smooth   = 0.23; % quality tolerance for smoothing interior nodes
SETTINGS.q_sliver   = 0.10; % quality tolerance for slivers (almost coplanar tetrahedrons)
SETTINGS.q_balloon  = 0.25; % quality tolerance for balloon forces
SETTINGS.mean_q_tol = 0.80; % quality tolerance for q mean
SETTINGS.mean_misfit_bar_length_tol = 0.14; % tolerance for the misfit of springs (it sets
                                            % when the right nodal density is achieved) 

% SOLVERS
SETTINGS.cross_bars            = 0; % whether to use cross bars (1--yes, 0--no)
                                    % (steadier structure in each triangle)
SETTINGS.cross_bars_faces      = 0; % whether to use cross bars on faces 
                                    % (1--yes, 0--no)
SETTINGS.balloon_forces        = 0; % whether to use ballon-forces (1--yes, 0--no)
                                    % (for having better triangles)
SETTINGS.solver_balloon_forces = 'sparse'; % solver used to compute incentre and inradius
                                           % 'gauss'     --> analytical solution
                                           %                 using Gauss Method 
                                           % '3D_matrix' --> 3-D matrix with "for"
                                           %                 loop 
                                           % 'sparse'    --> sparse matrix
SETTINGS.bc_method             = 'projection'; % keep projection!! uzawa not coded 
                                               % 'projection' --> fixing bnd nodes on tangent
                                               %                  line to circumference at
                                               %                  each node 
                                               % 'uzawa'      --> penalty method
SETTINGS.edges_output_mesh     = 'curved'; % set the shape of the edges for the quadratic 10nodel output mesh
                                           % 'curved'
                                           % 'straight'

%==========================================================================
% REFINEMENT SETTINGS
%==========================================================================
SETTINGS.refinement = 'guide_mesh'; % set refinement type
                                    % 'guide_mesh' --> refine the mesh using a guide-mesh
                                    % 'interface'  --> refine the mesh in a spherical layer
                                    % 'regular'    --> create a regular mesh
switch SETTINGS.refinement
    case 'guide_mesh'
        SETTINGS.guide_mesh = 'new'; % set either create or load a guide-mesh
                                     % 'new'  --> create a new guide mesh
                                     % 'load' --> load a pre-existing guide mesh
        SETTINGS.bnds_at_ref_region = 'fixed'; % set whether creating boundaries for the refined region.
                                               % 'fixed' -> Nodes confined in those boundaries are fixed
                                               % 'free'  -> Nodes confined in those boundaries will move 
                                               %            parallel to the boundaries 
                                               % 'no'    -> No boundaries for refined region
        switch SETTINGS.guide_mesh
            case 'new'
                % Point around which refined and transition zones are created
                GUIDE_MESH.theta0    = 90;     % colatitude (degrees) -> keep 90 for creating the mesh
                GUIDE_MESH.phi0      = 90;     % longitude (degrees)  -> keep 90 for creating the mesh
                GUIDE_MESH.r0        = 6371;   % radial distance (km)
                P_0                  = [GUIDE_MESH.theta0 GUIDE_MESH.phi0 GUIDE_MESH.r0];
                
                % Point around which the final mesh is centered (South Atlantic MOR 130 Myr ago) 
                GUIDE_MESH.theta_center = 128; % colatitude (degrees)
                GUIDE_MESH.phi_center   = 353; % longitude (degrees)
                P_center                = [GUIDE_MESH.theta_center GUIDE_MESH.phi_center GUIDE_MESH.r0];
                
                % Euler pole for the 2 finite rotation between (theta0,phi0) and (theta_center,phi_center).
                % We use 2 finite rotations in order to preserve the orientation of the refined zone 
                % 1st rotation around Z axis:
                P_mid      = P_center;
                P_mid(1,1) = P_0(1,1);
                [GUIDE_MESH.EP1_lat,GUIDE_MESH.EP1_lon,GUIDE_MESH.EP1_angle] = ...
                    euler_pole_from_two_points(P_0,P_mid);
                % 2nd rotation:
                [GUIDE_MESH.EP2_lat,GUIDE_MESH.EP2_lon,GUIDE_MESH.EP2_angle] = ...
                    euler_pole_from_two_points(P_mid,P_center);
                
                % % To compute the Euler Poles using GPlates:
                % %   1) Open GPlates
                % %   2) Click on Utilities -> Calculate Finite Rotation
                % %   3) Introduce initial and final points in the panel "Finite Rotation Between Points" 
                % GUIDE_MESH.EP1_lat   = -90;
                % GUIDE_MESH.EP1_lon   =   0;
                % GUIDE_MESH.EP1_angle =  97;
                % GUIDE_MESH.EP2_lat   =   0;
                % GUIDE_MESH.EP2_lon   =  83;
                % GUIDE_MESH.EP2_angle =  38;
                
                % Coarse zone
                GUIDE_MESH.l0_coarse = 2000;   % desired spring length (km) for coarse zone
                
                % Transition zone
                GUIDE_MESH.d_tran    = 2900; % transition zone depth (km)
                GUIDE_MESH.w_tran    = 9600; %8000; % transition zone width (km) (North-South)
                GUIDE_MESH.l_tran    = 6800; %5600; % transition zone length (km) (East-West)
                
                % Refined zone
                GUIDE_MESH.d_ref     = 300;  % refined zone depth (km)
                GUIDE_MESH.w_ref     = 5000; %4200; % refined zone width (km) (North-South)
                GUIDE_MESH.l_ref     = 2200; %1800; % refined zone length (km) (East-West)
                GUIDE_MESH.l0_ref    = 60;   %100;  % desired spring length (km) for refined zone
                
            case 'load'
                % Full path to file containing the guide mesh is required!
                SETTINGS.guide_mesh_file = 'SETUP_TEST/GUIDE_MESH';
                GUIDE_MESH               = struct([]); % empty structure
            otherwise
                error('typo in SETTINGS.guide_mesh')
        end
        INTERFACE = struct([]); % empty structure
    
    case 'interface'
        SETTINGS.h0          = 3.5; % desired spring length (km) for a regular mesh
        INTERFACE.r          = 5;   % radius of the interface 
        INTERFACE.ref_factor = 16;  % factor between the regular elements and refined ones
        INTERFACE.s          = 0.8; % slope to smoothly increase the element size
        GUIDE_MESH           = struct([]); % empty structure
    
    case 'regular'
        SETTINGS.h0          = 1000;       % desired spring length (km) for a regular mesh
        GUIDE_MESH           = struct([]); % empty structure
        INTERFACE            = struct([]); % empty structure
    
    otherwise
        error('typo in SETTINGS.refinement')
end

%==========================================================================
% FIRST GUESS SETTINGS
%==========================================================================
SETTINGS.first_guess = 'load';
    % 'new'          --> create a new first guess for the mesh
    % 'load'         --> load a pre-existing mesh
    % 'load_as_pfix' --> load a pre-existing mesh as pfix and create the rest of nodes 
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
path2data       = data_storage_3d(); % edit file "data_storage_3d" to add 
                                     % location on new computer
SETTINGS.outdir = [path2data 'SPH_MESH/Trash_01']; % full path to output 

end % END OF FUNCTION mesh_parameters