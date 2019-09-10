function [pfix,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH)
% Usage: [pfix,pbnd1,pbnd2] = boundary_nodes(SETTINGS,GUIDE_MESH)
%
% Purpose:
%   Generation of boundary nodes for a 2D cylindrical annulus mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   pfix       : [matrix]    : coord of fixed nodes
%   pbnd1      : [matrix]    : coord of boundary nodes for inner boundary
%   pbnd2      : [matrix]    : coord of boundary nodes for outer boundary
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
deg2rad      = pi/180;
r_int        = SETTINGS.r_int;      % inner radius (km) of cylindrical annulus
r_ext        = SETTINGS.r_ext;      % outer radius (km) of cylindrical annulus
theta0       = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
w_tran_deg   = GUIDE_MESH.w_tran/(deg2rad*r_ext); % width of transition zone in degrees
theta_tran_l = theta0-w_tran_deg/2; % colatitude of the left boundary in the transition zone
theta_tran_r = theta0+w_tran_deg/2; % colatitude of the right boundary in the transition zone
w_ref_deg    = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees
theta_ref_l  = theta0-w_ref_deg/2;  % colatitude of the left boundary in the refined zone
theta_ref_r  = theta0+w_ref_deg/2;  % colatitude of the right boundary in the refined zone

%==========================================================================
% BOUNDARY 1 (INNER BOUNDARY)
%==========================================================================
% COARSE ZONE
pbnd1_coarse         = regular_bnd_nodes(r_int,GUIDE_MESH.l0_coarse); 
[pfix1,pbnd1_coarse] = create_pfix(pbnd1_coarse,r_int);
[pbnd1_coarse_pol]   = cartesian2polar(pbnd1_coarse);
pbnd1_coarse(pbnd1_coarse_pol(:,1) > theta_tran_l &...
             pbnd1_coarse_pol(:,1) < theta_tran_r,:) = []; % remove bnd nodes that are in the transition zone
pbnd1_coarse         = reject_points_line(pbnd1_coarse,GUIDE_MESH); % reject pbnd1_coarse using probability

% TRANSITION ZONE
pbnd1_tran           = regular_bnd_nodes(r_int,GUIDE_MESH.l0_coarse);
[~,pbnd1_tran]       = create_pfix(pbnd1_tran,r_int); % remove pfix1 values from pbnd1_tran
pbnd1_tran_pol       = cartesian2polar(pbnd1_tran);
theta1_tran          = pbnd1_tran_pol(:,1);
pbnd1_tran           = pbnd1_tran(theta1_tran > theta_tran_l &...
                                  theta1_tran < theta_tran_r,:); % take only the nodes in the transition zone
pbnd1_tran           = reject_points_line(pbnd1_tran,GUIDE_MESH); % reject pbnd1_tran using probability

%==========================================================================
% BOUNDARY 2 (OUTER BOUNDARY)
%==========================================================================
% COARSE ZONE
pbnd2_coarse         = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_coarse);
[~,pbnd2_coarse] = create_pfix(pbnd2_coarse,r_ext);
[pbnd2_coarse_pol]   = cartesian2polar(pbnd2_coarse);
pbnd2_coarse(pbnd2_coarse_pol(:,1) > theta_tran_l &...
             pbnd2_coarse_pol(:,1) < theta_tran_r,:) = []; % remove bnd nodes that are in the transition zone
pbnd2_coarse         = reject_points_line(pbnd2_coarse,GUIDE_MESH); % reject pbnd2_coarse using probability

% TRANSITION ZONE
pbnd2_tran           = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_ref);
[~,pbnd2_tran]       = create_pfix(pbnd2_tran,r_ext); % remove pfix2 values from pbnd2_tran
pbnd2_tran_pol       = cartesian2polar(pbnd2_tran);
theta2_tran          = pbnd2_tran_pol(:,1);
pbnd2_tran           = pbnd2_tran(theta2_tran > theta_tran_l &...
                                  theta2_tran < theta_tran_r,:); % take only the nodes in the transition zone
theta2_tran_temp     = theta2_tran(theta2_tran > theta_tran_l &...
                                   theta2_tran < theta_tran_r,:); % take theta angles inside the transition zone
pbnd2_tran(theta2_tran_temp > theta_ref_l &...
           theta2_tran_temp < theta_ref_r,:) = []; % remove bnd nodes that are inside the refined zone
pbnd2_tran           = reject_points_line(pbnd2_tran,GUIDE_MESH); % reject pbnd2_tran using probability
       
% REFINED ZONE
pbnd2_ref            = regular_bnd_nodes(r_ext,GUIDE_MESH.l0_ref);
[pfix2,pbnd2_ref]    = create_pfix(pbnd2_ref,r_ext); % remove pfix2 values from pbnd2_ref
pbnd2_ref_pol        = cartesian2polar(pbnd2_ref);
theta2_ref           = pbnd2_ref_pol(:,1);
pbnd2_ref            = pbnd2_ref(theta2_ref > theta_ref_l &...
                                 theta2_ref < theta_ref_r,:); % take only the nodes in the refined zone
pbnd2_ref            = reject_points_line(pbnd2_ref,GUIDE_MESH); % reject pbnd2_ref using probability

%==========================================================================
% OUPUT DATA
%==========================================================================
pfix  = [pfix1; pfix2];
pbnd1 = [pbnd1_coarse; pbnd1_tran];
pbnd2 = [pbnd2_coarse; pbnd2_tran; pbnd2_ref];

if strcmp(SETTINGS.mesh,'axisym')
    if SETTINGS.r_int == 0;
        y_ref      = linspace(GUIDE_MESH.r0, GUIDE_MESH.r0 - GUIDE_MESH.d_ref, round(GUIDE_MESH.d_ref/GUIDE_MESH.l0_ref))';
        y_tran     = [GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 30   ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 60   ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 100  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 150  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 200  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 260  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 320  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 390  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 460  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 540  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 620  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 720  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 840  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 980  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1140 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1320 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1540 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1810 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 2150 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 2600];
        y_coarse   = linspace(GUIDE_MESH.r0 - GUIDE_MESH.d_tran - GUIDE_MESH.l0_coarse, -r_ext, ...
                              round((GUIDE_MESH.r0 - GUIDE_MESH.d_tran - GUIDE_MESH.l0_coarse + r_ext)/GUIDE_MESH.l0_coarse))';
        y          = [y_ref; y_tran; y_coarse];
        y([1,end]) = [];
        x          = zeros(size(y,1),1);
        pfix(3,:)  = [0 -r_ext];
        pfix       = [pfix; [x y]];
        pfix(pfix(:,1) < 0,:)    = []; % remove pfix points with x < 0
        pbnd2(pbnd2(:,1) <= 0,:) = []; % remove pbnd2 points with x <= 0
    else
        y_ref      = linspace(GUIDE_MESH.r0,GUIDE_MESH.r0-GUIDE_MESH.d_ref,round(GUIDE_MESH.d_ref/GUIDE_MESH.l0_ref))';
        y_tran     = [GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 30   ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 70   ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 125  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 195  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 285  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 395  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 530  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 705  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 905  ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1135 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1400 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 1700 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 2100 ; ...
                      GUIDE_MESH.r0 - GUIDE_MESH.d_ref - 2600];
        y1          = [y_ref; y_tran];
        y1([1,end]) = [];
        x1          = zeros(size(y1,1),1);
        y2          = linspace(-r_int,-r_ext,round((r_ext-r_int)/GUIDE_MESH.l0_coarse))';
        y2([1,end]) = [];
        x2          = zeros(size(y2,1),1);
        pfix(3,:)   = [0 -r_int];
        pfix(7,:)   = [0 -r_ext];
        pfix        = [pfix; [x1 y1]; [x2 y2]];
        pfix(pfix(:,1) < 0,:)    = []; % remove pfix points with x < 0
        pbnd1(pbnd1(:,1) <= 0,:) = []; % remove pbnd1 points with x <= 0
        pbnd2(pbnd2(:,1) <= 0,:) = []; % remove pbnd2 points with x <= 0
    end
end

end % END OF FUNCTION boundary_nodes

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function GCOORD = reject_points_line(GCOORD,GUIDE_MESH)
% Usage: GCOORD = reject_points_line(GCOORD,GUIDE_MESH)
%
% Purpose:
%   Reject points on a line using point density (based on Persson and
%   Strang, 2004)
%
% Input:
%   GCOORD     : [matrix]    : coordinates of points
%   GUIDE_MESH : [structure] : structure containing guide mesh settings
%
% Output:
%   GCOORD     : [matrix]    : coordinates of points
%
% JMT Jan 2017

L0_GCOORD = bar_L0_guide(GCOORD,GUIDE_MESH); % compute the desired length for each point
r0_GCOORD = sqrt(2)./L0_GCOORD;              % probability to keep the points (it is proportional to L^-1 because points are on a line)
GCOORD    = GCOORD(rand(size(GCOORD,1),1)<r0_GCOORD./max(r0_GCOORD),:); % reject points with a probability proportional to L^-1

end % END OF FUNCTION reject_points_line