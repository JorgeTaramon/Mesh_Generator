function MESH = projection_v2(MESH,GUIDE_MESH,STIFFNESS,F_le,SETTINGS)
% Usage: MESH = projection_v2(MESH,GUIDE_MESH,STIFFNESS,F_le,SETTINGS)
%
% Purpose:
%   Fix boundary nodes motion on the tangent plane to each boundary node 
%   and solve for the answer.
%   We rotate bnd nodes to their own axis domain (X",Y",Z"), while the 
%   tangent plane is X"-Y". So, by fixing Z" direction displacements, we 
%   fix the bnd nodes on their tangent planes to the sphere.
%   We design transformation matrix [TT] for this goal which is in function
%   of phi (longitude, measured from X axis in counterclowise direction), 
%   and theta (colatitude, measured from +Z axis to -Z axis).
%   The transformation matrix is the result after applying one rotation 
%   around the Z axis and another one around the Y' axis.
%
%                 Z
%                 |
%                 |_
%                 | \
%                 |  \
%                 |   \
%                 |    \
%                 |     +
%                 |    /.
%                 |th / .
%                 |  /r .
%                 |^/   .
%                 |/    .
%                _-_----.-------_----> Y
%              _-\_/\   .     _-
%            _-  phi \  .   _-
%          _-         \ . _-
%        _-____________\.-
%      _-                
%    X      
%
%   1st rotation (Counterlockwise rotation around Z axis an angle phi):
%
%   [x ]   [ cos(phi)  -sin(phi)   0  ][x']
%   |y | = | sin(phi)   cos(phi)   0  ||y'|
%   [z ]   [    0          0       1  ][z']
%
%   2nd rotation (Counterlockwise around Y' axis an angle theta):
%
%   [x']   [ cos(theta)  0  sin(theta)][x"]
%   |y'| = |     0       1      0     ||y"|
%   [z']   [-sin(theta)  0  cos(theta)][z"]
%
%   The transformation matrix [TT] to change from local coordinates
%   (x",y",z") to global coordinates (x,y,z) is given by:
%   
%   [x]   [cos(phi)cos(theta)  -sin(phi)   cos(phi)sin(theta)][x"]
%   |y| = |sin(phi)cos(theta)   cos(phi)   sin(phi)sin(theta)||y"|
%   [z]   [   -sin(theta)         0            cos(theta)    ][z"] (= r for bnd nodes) <-- constrain BC
%
% Input:
%   MESH       : [structure]     : structure containing the mesh
%   GUIDE_MESH : [structure]     : structure containing guide mesh
%   STIFFNESS  : [structure]     : structure containing matrix parameters
%   F_le       : [column vector] : force vector term due to the effect
%                                  of a preferred zero-deformation element
%                                  length
%   SETTINGS   : [structure]     : structure containing mesh settings
%
% Output:
%   MESH       : [structure]     : structure containing the mesh
%
% JMT Jul 2015
% JMT May 2016: cleaned up
% JMT Aug 2017: Handles internal boundaries. The nodes of the internal
%               boundaries for the refined region are free to move along 
%               the internal boundaries.
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
r_int             = SETTINGS.r_int;
r_ext             = SETTINGS.r_ext;
GCOORD            = MESH.GCOORD;
nnod              = size(GCOORD,1);
nfix              = size(MESH.pfix,1);
nbnd1             = size(MESH.pbnd1,1);
nbnd2             = size(MESH.pbnd2,1);
nbnd_ref_face_bot = size(MESH.pbnd_ref_face_bot,1);
nbnd_ref_face_1   = size(MESH.pbnd_ref_face_1,1);
nbnd_ref_face_2   = size(MESH.pbnd_ref_face_2,1);
nbnd_ref_face_3   = size(MESH.pbnd_ref_face_3,1);
nbnd_ref_face_4   = size(MESH.pbnd_ref_face_4,1);
GCOORD_init       = GCOORD;           % save for plotting
deg2rad           = pi/180;
theta0            = GUIDE_MESH.theta0;   % colatitude (degrees) of the point around which the refined and transition zones are defined
phi0              = GUIDE_MESH.phi0;     % longitude (degrees) of the point around which the refined and transition zones are defined
d_ref             = GUIDE_MESH.d_ref;    % refined zone depth (km)
w_ref_deg         = GUIDE_MESH.w_ref/(deg2rad*r_ext);  % width of refined zone in degrees (North-South)
theta_ref_n       = theta0-w_ref_deg/2;  % colatitude of the northern boundary in the refined zone
theta_ref_s       = theta0+w_ref_deg/2;  % colatitude of the southern boundary in the refined zone
l_ref_deg         = GUIDE_MESH.l_ref/(deg2rad*r_ext);  % length of refined zone in degrees (East-West)
phi_ref_e         = phi0+l_ref_deg/2;    % longitude of the eastern boundary in the refined zone
phi_ref_w         = phi0-l_ref_deg/2;    % longitude of the western boundary in the refined zone

%==========================================================================
% ASSEMBLE TRANSFORMATION MATRIX [TT]
%==========================================================================
% theta and phi angles for nodes on both boundaries
if r_int > 0 && ~isempty(MESH.pbnd1)
    pbnd1_sph      = cartesian2spherical(MESH.pbnd1);
    theta1         = pbnd1_sph(:,1);
    phi1           = pbnd1_sph(:,2);
end
pbnd2_sph       = cartesian2spherical(MESH.pbnd2);
theta2          = pbnd2_sph(:,1);
phi2            = pbnd2_sph(:,2);

% theta and phi angles for nodes on internal faces
pbnd_ref_face_bot_sph = cartesian2spherical(MESH.pbnd_ref_face_bot);
th_face_bot           = pbnd_ref_face_bot_sph(:,1);
ph_face_bot           = pbnd_ref_face_bot_sph(:,2);
pbnd_ref_face_1_sph   = cartesian2spherical(MESH.pbnd_ref_face_1);
th_face_1             = pbnd_ref_face_1_sph(:,1);
ph_face_1             = pbnd_ref_face_1_sph(:,2);
pbnd_ref_face_2_sph   = cartesian2spherical(MESH.pbnd_ref_face_2);
th_face_2             = pbnd_ref_face_2_sph(:,1);
ph_face_2             = pbnd_ref_face_2_sph(:,2);
pbnd_ref_face_3_sph   = cartesian2spherical(MESH.pbnd_ref_face_3);
th_face_3             = pbnd_ref_face_3_sph(:,1);
ph_face_3             = pbnd_ref_face_3_sph(:,2);
pbnd_ref_face_4_sph   = cartesian2spherical(MESH.pbnd_ref_face_4);
th_face_4             = pbnd_ref_face_4_sph(:,1);
ph_face_4             = pbnd_ref_face_4_sph(:,2);

% pre-allocate memory: number of non-zero slots in [TT]
nnz_T           = 3*nnod + 6*nbnd1 + 6*nbnd2 + ...
                  6*nbnd_ref_face_bot + ...
                  6*nbnd_ref_face_1 + 6*nbnd_ref_face_2 + ...
                  6*nbnd_ref_face_3 + 6*nbnd_ref_face_4; 
                % = 3*nfix + 9*nbnd1 + 9*nbnd2 +            ...
                %   9*nbnd_ref_face_bot +                   ...
                %   9*nbnd_ref_face_1 + 9*nbnd_ref_face_1 + ...
                %   9*nbnd_ref_face_1 + 9*nbnd_ref_face_1 + 3*nint
KKi_T           = zeros(nnz_T,1); 
KKj_T           = KKi_T; 
KK1_T           = ones(nnz_T,1);
KKi_T(1:3*nnod) = (1:3*nnod)';
KKj_T(1:3*nnod) = (1:3*nnod)';

if r_int > 0 && ~isempty(MESH.pbnd1)
    % diagonal slots for boundary 1: cos(phi1)*cos(theta1), cos(phi1), cos(theta1) for each node in bnd 1
    % (3 slots for each bnd node).
    KK1_T(3*nfix+1:3:3*(nfix+nbnd1)-2) = cosd(phi1).*cosd(theta1); % element (1,1) in [TT] for bnd 1 nodes
    KK1_T(3*nfix+2:3:3*(nfix+nbnd1)-1) = cosd(phi1);               % element (2,2) in [TT] for bnd 1 nodes
    KK1_T(3*nfix+3:3:3*(nfix+nbnd1))   = cosd(theta1);             % element (3,3) in [TT] for bnd 1 nodes
end

% diagonal slots for boundary 2: cos(phi2)*cos(theta2), cos(phi2), cos(theta2) for each node in bnd 2
% (3 slots for each bnd node).
KK1_T(3*(nfix+nbnd1)+1:3:3*(nfix+nbnd1+nbnd2)-2) = cosd(phi2).*cosd(theta2); % element (1,1) in [TT] for bnd 2 nodes
KK1_T(3*(nfix+nbnd1)+2:3:3*(nfix+nbnd1+nbnd2)-1) = cosd(phi2);               % element (2,2) in [TT] for bnd 2 nodes
KK1_T(3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2))   = cosd(theta2);             % element (3,3) in [TT] for bnd 2 nodes

% diagonal slots for pint_ref_face_bot:
KK1_T(3*(nfix+nbnd1+nbnd2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2) = cosd(ph_face_bot).*cosd(th_face_bot); % element (1,1) in [TT] for pbnd_ref_face_bot
KK1_T(3*(nfix+nbnd1+nbnd2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1) = cosd(ph_face_bot);                    % element (2,2) in [TT] for pbnd_ref_face_bot
KK1_T(3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot))   = cosd(th_face_bot);                    % element (3,3) in [TT] for pbnd_ref_face_bot

% diagonal slots for pbnd_ref_face_1:
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2) = cosd(ph_face_1).*cosd(th_face_1); % element (1,1) in [TT] for pbnd_ref_face_1
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1) = cosd(ph_face_1);                  % element (2,2) in [TT] for pbnd_ref_face_1
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1))   = cosd(th_face_1);                  % element (3,3) in [TT] for pbnd_ref_face_1

% diagonal slots for pbnd_ref_face_2:
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2) = cosd(ph_face_2).*cosd(th_face_2); % element (1,1) in [TT] for pbnd_ref_face_2
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1) = cosd(ph_face_2);                  % element (2,2) in [TT] for pbnd_ref_face_2
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2))   = cosd(th_face_2);                  % element (3,3) in [TT] for pbnd_ref_face_2

% diagonal slots for pbnd_ref_face_3:
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2) = cosd(ph_face_3).*cosd(th_face_3); % element (1,1) in [TT] for pbnd_ref_face_3
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1) = cosd(ph_face_3);                  % element (2,2) in [TT] for pbnd_ref_face_3
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3))   = cosd(th_face_3);                  % element (3,3) in [TT] for pbnd_ref_face_3

% diagonal slots for pbnd_ref_face_4:
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2) = cosd(ph_face_4).*cosd(th_face_4); % element (1,1) in [TT] for pbnd_ref_face_4
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1) = cosd(ph_face_4);                  % element (2,2) in [TT] for pbnd_ref_face_4
KK1_T(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4))   = cosd(th_face_4);                  % element (3,3) in [TT] for pbnd_ref_face_4

if r_int > 0 && ~isempty(MESH.pbnd1)
    % off-diagonal slots for boundary 1. This will use up 6*nbnd1 slots
    KKi_T(3*nnod+        1:3*nnod+  nbnd1) = (3*nfix+1:3:3*(nfix+nbnd1)-2);
    KKj_T(3*nnod+        1:3*nnod+  nbnd1) = (3*nfix+2:3:3*(nfix+nbnd1)-1);
    KK1_T(3*nnod+        1:3*nnod+  nbnd1) = -sind(phi1);              % element (1,2) in [TT] for bnd 1 nodes
    KKi_T(3*nnod+  nbnd1+1:3*nnod+2*nbnd1) = (3*nfix+1:3:3*(nfix+nbnd1)-2);
    KKj_T(3*nnod+  nbnd1+1:3*nnod+2*nbnd1) = (3*nfix+3:3:3*(nfix+nbnd1)  );
    KK1_T(3*nnod+  nbnd1+1:3*nnod+2*nbnd1) = cosd(phi1).*sind(theta1); % element (1,3) in [TT] for bnd 1 nodes
    KKi_T(3*nnod+2*nbnd1+1:3*nnod+3*nbnd1) = (3*nfix+2:3:3*(nfix+nbnd1)-1);
    KKj_T(3*nnod+2*nbnd1+1:3*nnod+3*nbnd1) = (3*nfix+1:3:3*(nfix+nbnd1)-2);
    KK1_T(3*nnod+2*nbnd1+1:3*nnod+3*nbnd1) = sind(phi1).*cosd(theta1); % element (2,1) in [TT] for bnd 1 nodes
    KKi_T(3*nnod+3*nbnd1+1:3*nnod+4*nbnd1) = (3*nfix+2:3:3*(nfix+nbnd1)-1);
    KKj_T(3*nnod+3*nbnd1+1:3*nnod+4*nbnd1) = (3*nfix+3:3:3*(nfix+nbnd1)  );
    KK1_T(3*nnod+3*nbnd1+1:3*nnod+4*nbnd1) = sind(phi1).*sind(theta1); % element (2,3) in [TT] for bnd 1 nodes
    KKi_T(3*nnod+4*nbnd1+1:3*nnod+5*nbnd1) = (3*nfix+3:3:3*(nfix+nbnd1)  );
    KKj_T(3*nnod+4*nbnd1+1:3*nnod+5*nbnd1) = (3*nfix+1:3:3*(nfix+nbnd1)-2);
    KK1_T(3*nnod+4*nbnd1+1:3*nnod+5*nbnd1) = -sind(theta1);            % element (3,1) in [TT] for bnd 1 nodes
    KKi_T(3*nnod+5*nbnd1+1:3*nnod+6*nbnd1) = (3*nfix+3:3:3*(nfix+nbnd1)  );
    KKj_T(3*nnod+5*nbnd1+1:3*nnod+6*nbnd1) = (3*nfix+2:3:3*(nfix+nbnd1)-1);
    KK1_T(3*nnod+5*nbnd1+1:3*nnod+6*nbnd1) = 0;                        % element (3,2) in [TT] for bnd 1 nodes
end

% off-diagonal slots for boundary 2. This will use up 6*nbnd2 slots
KKi_T(3*nnod+6*nbnd1+        1:3*nnod+6*nbnd1+  nbnd2) = (3*(nfix+nbnd1)+1:3:3*(nfix+nbnd1+nbnd2)-2);
KKj_T(3*nnod+6*nbnd1+        1:3*nnod+6*nbnd1+  nbnd2) = (3*(nfix+nbnd1)+2:3:3*(nfix+nbnd1+nbnd2)-1);
KK1_T(3*nnod+6*nbnd1+        1:3*nnod+6*nbnd1+  nbnd2) = -sind(phi2);              % element (1,2) in [TT] for bnd 2 nodes
KKi_T(3*nnod+6*nbnd1+  nbnd2+1:3*nnod+6*nbnd1+2*nbnd2) = (3*(nfix+nbnd1)+1:3:3*(nfix+nbnd1+nbnd2)-2);
KKj_T(3*nnod+6*nbnd1+  nbnd2+1:3*nnod+6*nbnd1+2*nbnd2) = (3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2)  );
KK1_T(3*nnod+6*nbnd1+  nbnd2+1:3*nnod+6*nbnd1+2*nbnd2) = cosd(phi2).*sind(theta2); % element (1,3) in [TT] for bnd 2 nodes
KKi_T(3*nnod+6*nbnd1+2*nbnd2+1:3*nnod+6*nbnd1+3*nbnd2) = (3*(nfix+nbnd1)+2:3:3*(nfix+nbnd1+nbnd2)-1);
KKj_T(3*nnod+6*nbnd1+2*nbnd2+1:3*nnod+6*nbnd1+3*nbnd2) = (3*(nfix+nbnd1)+1:3:3*(nfix+nbnd1+nbnd2)-2);
KK1_T(3*nnod+6*nbnd1+2*nbnd2+1:3*nnod+6*nbnd1+3*nbnd2) = sind(phi2).*cosd(theta2); % element (2,1) in [TT] for bnd 2 nodes
KKi_T(3*nnod+6*nbnd1+3*nbnd2+1:3*nnod+6*nbnd1+4*nbnd2) = (3*(nfix+nbnd1)+2:3:3*(nfix+nbnd1+nbnd2)-1);
KKj_T(3*nnod+6*nbnd1+3*nbnd2+1:3*nnod+6*nbnd1+4*nbnd2) = (3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2)  );
KK1_T(3*nnod+6*nbnd1+3*nbnd2+1:3*nnod+6*nbnd1+4*nbnd2) = sind(phi2).*sind(theta2); % element (2,3) in [TT] for bnd 2 nodes
KKi_T(3*nnod+6*nbnd1+4*nbnd2+1:3*nnod+6*nbnd1+5*nbnd2) = (3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2)  );
KKj_T(3*nnod+6*nbnd1+4*nbnd2+1:3*nnod+6*nbnd1+5*nbnd2) = (3*(nfix+nbnd1)+1:3:3*(nfix+nbnd1+nbnd2)-2);
KK1_T(3*nnod+6*nbnd1+4*nbnd2+1:3*nnod+6*nbnd1+5*nbnd2) = -sind(theta2);            % element (3,1) in [TT] for bnd 2 nodes
KKi_T(3*nnod+6*nbnd1+5*nbnd2+1:3*nnod+6*nbnd1+6*nbnd2) = (3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2)  );
KKj_T(3*nnod+6*nbnd1+5*nbnd2+1:3*nnod+6*nbnd1+6*nbnd2) = (3*(nfix+nbnd1)+2:3:3*(nfix+nbnd1+nbnd2)-1);
KK1_T(3*nnod+6*nbnd1+5*nbnd2+1:3*nnod+6*nbnd1+6*nbnd2) = 0;                        % element (3,2) in [TT] for bnd 2 nodes

% off-diagonal slots for pbnd_ref_face_bot. This will use up 6*nbnd_ref_face_bot slots
KKi_T(3*nnod+6*nbnd1+6*nbnd2+                    1:3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+                    1:3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+                    1:3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot) = -sind(ph_face_bot);                   % element (1,2) in [TT] for pbnd_ref_face_bot
KKi_T(3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+  nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot) = cosd(ph_face_bot).*sind(th_face_bot); % element (1,3) in [TT] for pbnd_ref_face_bot
KKi_T(3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+2*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot) = sind(ph_face_bot).*cosd(th_face_bot); % element (2,1) in [TT] for pbnd_ref_face_bot
KKi_T(3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+3*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot) = sind(ph_face_bot).*sind(th_face_bot); % element (2,3) in [TT] for pbnd_ref_face_bot
KKi_T(3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+4*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot) = -sind(th_face_bot);                   % element (3,1) in [TT] for pbnd_ref_face_bot
KKi_T(3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot) = (3*(nfix+nbnd1+nbnd2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+5*nbnd_ref_face_bot+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot) = 0;                                    % element (3,2) in [TT] for pbnd_ref_face_bot

% off-diagonal slots for pbnd_ref_face_1. This will use up 6*nbnd_ref_face_1 slots
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1) = -sind(ph_face_1);                 % element (1,2) in [TT] for pbnd_ref_face_1
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+  nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1) = cosd(ph_face_1).*sind(th_face_1); % element (1,3) in [TT] for pbnd_ref_face_1
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+2*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1) = sind(ph_face_1).*cosd(th_face_1); % element (2,1) in [TT] for pbnd_ref_face_1
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+3*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1) = sind(ph_face_1).*sind(th_face_1); % element (2,3) in [TT] for pbnd_ref_face_1
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+4*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1) = -sind(th_face_1);                 % element (3,1) in [TT] for pbnd_ref_face_1
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+5*nbnd_ref_face_1+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1) = 0;                                % element (3,2) in [TT] for pbnd_ref_face_1

% off-diagonal slots for pbnd_ref_face_2. This will use up 6*nbnd_ref_face_2 slots
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2) = -sind(ph_face_2);                 % element (1,2) in [TT] for pbnd_ref_face_2
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+  nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2) = cosd(ph_face_2).*sind(th_face_2); % element (1,3) in [TT] for pbnd_ref_face_2
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+2*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2) = sind(ph_face_2).*cosd(th_face_2); % element (2,1) in [TT] for pbnd_ref_face_2
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+3*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2) = sind(ph_face_2).*sind(th_face_2); % element (2,3) in [TT] for pbnd_ref_face_2
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+4*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2) = -sind(th_face_2);                 % element (3,1) in [TT] for pbnd_ref_face_2
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+5*nbnd_ref_face_2+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2) = 0;                                % element (3,2) in [TT] for pbnd_ref_face_2

% off-diagonal slots for pbnd_ref_face_3. This will use up 6*nbnd_ref_face_3 slots
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3) = -sind(ph_face_3);                 % element (1,2) in [TT] for pbnd_ref_face_3
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+  nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3) = cosd(ph_face_3).*sind(th_face_3); % element (1,3) in [TT] for pbnd_ref_face_3
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+2*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3) = sind(ph_face_3).*cosd(th_face_3); % element (2,1) in [TT] for pbnd_ref_face_3
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+3*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3) = sind(ph_face_3).*sind(th_face_3); % element (2,3) in [TT] for pbnd_ref_face_3
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+4*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3) = -sind(th_face_3);                 % element (3,1) in [TT] for pbnd_ref_face_3
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+5*nbnd_ref_face_3+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3) = 0;                                % element (3,2) in [TT] for pbnd_ref_face_3

% off-diagonal slots for pbnd_ref_face_4. This will use up 6*nbnd_ref_face_4 slots
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+                  1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4) = -sind(ph_face_4);                 % element (1,2) in [TT] for pbnd_ref_face_4
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+  nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4) = cosd(ph_face_4).*sind(th_face_4); % element (1,3) in [TT] for pbnd_ref_face_4
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+2*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4) = sind(ph_face_4).*cosd(th_face_4); % element (2,1) in [TT] for pbnd_ref_face_4
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1);
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)  );
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+3*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4) = sind(ph_face_4).*sind(th_face_4); % element (2,3) in [TT] for pbnd_ref_face_4
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4) = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+4*nbnd_ref_face_4+1:3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4) = -sind(th_face_4);                 % element (3,1) in [TT] for pbnd_ref_face_4
KKi_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4+1:end)                                                                                                                = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)  );
KKj_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4+1:end)                                                                                                                = (3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1);
KK1_T(3*nnod+6*nbnd1+6*nbnd2+6*nbnd_ref_face_bot+6*nbnd_ref_face_1+6*nbnd_ref_face_2+6*nbnd_ref_face_3+5*nbnd_ref_face_4+1:end)                                                                                                                = 0;                                % element (3,2) in [TT] for pbnd_ref_face_4

TT = sparse2(KKi_T(:),KKj_T(:),KK1_T(:)); % transformation matrix

%==========================================================================
% ASSEMBLE THE GLOBAL STIFFNESS MATRIX [KK] AND SOLVE FOR ANSWER {X}
%==========================================================================
% The system of equations to be solved before applying the spherical BCs is
% given by:
%
% [KK]{X} = {F} + {F_le} where:
%   {F} = 0 since we are solving for steady state of the springs
%   {F_le} is the force vector term due to the effect of a preferred
%    zero-deformation element length
%
% After applying the transformation matrix [TT] (to impose the spherical
% BCs) to each term
%
% [KK"]   = [TT]'[KK][TT]
% {X"}    = [TT]'{X}
% {F_le"} = [TT]'{F_le}
%
% the system of equations became:
%
% [KK"]{X"} = {F_le"}
%
% which is solved for {X"} where the BCs have been imposed and passed to 
% the right hand side as constants (fix).
% Finally, the global coordinates {X} can be recovered through the
% transformation matrix [TT]:
%
% {X} = [TT]{X"}

KKi = [STIFFNESS.KKi(:);STIFFNESS.KKi_CrossBars(:);STIFFNESS.KKi_CrossBars_faces(:)];
KKj = [STIFFNESS.KKj(:);STIFFNESS.KKj_CrossBars(:);STIFFNESS.KKj_CrossBars_faces(:)];
KK1 = [STIFFNESS.KK1(:);STIFFNESS.KK1_CrossBars(:);STIFFNESS.KK1_CrossBars_faces(:)];

KK_0 = sparse2(KKi(:),KKj(:),KK1(:)); % MATLAB's sparse matrix maker accumulates all element KK's
free          = [3*nfix+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-2, ... % theta is free (bnd1,bnd2 and bnd_ref_face_bot)
                 3*nfix+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)-1, ... % phi is free (bnd1,bnd2 and bnd_ref_face_bot)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-1, ... % phi is free (bnd_ref_face_1)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)  , ... % r is free (bnd_ref_face_1)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-2, ... % theta is free (bnd_ref_face_2)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)  , ... % r is free (bnd_ref_face_2)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-1, ... % phi is free (bnd_ref_face_3)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)  , ... % r is free (bnd_ref_face_3)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-2, ... % theta is free (bnd_ref_face_4)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)  , ... % r is free (bnd_ref_face_4)
                 3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)+1:3*nnod]; % theta, phi and r are free indices (interior nodes)
all_ind       = 1:3*nnod;
all_ind(free) = [];
fix           = all_ind; % fix indices

% Apply transformation matrix (to impose spherical BCs) to stiffness
% matrix [KK"] = [TT]'[KK][TT] in two steps:
% 1st step: [KK_1] = [KK_0][TT]
KK_1  = KK_0*TT;


% Apply transformation matrix (transposed) to global coordinates {X"}=[TT]'{X}
% (to impose spherical BCs in local coordinates)
GCOORD_temp = GCOORD';
GCOORD_temp = GCOORD_temp(:);
GCOORD_temp = TT'*GCOORD_temp;
if r_int > 0 && ~isempty(MESH.pbnd1)
    % Impose BCs for boundary 1, fixing the z" coordinate (of bnd nodes) with r_int
    GCOORD_temp(3*nfix+3:3:3*(nfix+nbnd1)) = r_int;
end
% Impose BCs for boundary 2, fixing the z" coordinate (of bnd nodes) with r_ext
GCOORD_temp(3*(nfix+nbnd1)+3:3:3*(nfix+nbnd1+nbnd2)) = r_ext;
% Impose BCs for pbnd_ref_face_bot, fixing the z" coordinate with r_ext - d_ref 
GCOORD_temp(3*(nfix+nbnd1+nbnd2)+3:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)) = r_ext - d_ref;
% Impose BCs for pbnd_ref_face_1, fixing the x" coordinate with 0 
GCOORD_temp(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)-2) = 0;
% Impose BCs for pbnd_ref_face_2, fixing the y" coordinate with 0
GCOORD_temp(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)-1) = 0;
% Impose BCs for pbnd_ref_face_3, fixing the x" coordinate with 0 
GCOORD_temp(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2)+1:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)-2) = 0;
% Impose BCs for pbnd_ref_face_4, fixing the y" coordinate with 0
GCOORD_temp(3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3)+2:3:3*(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4)-1) = 0;

% Pass imposed BCs to the right hand side in two steps:
% 1st step: rhs = -([KK_0][TT])[TT]'{X} + {F_le}
rhs = -KK_1(:,fix)*GCOORD_temp(fix) + F_le ; % = 0 - KK*pfix
% 2nd step: rhs = [TT]'(-([KK_0][TT])[TT]'{X} + {F_le}) =
%
%               = -([TT]'[KK_0][TT])[TT]'{X} + [TT]'{F_le} =
%
%               = - [KK"(fix)]{X"(fix)} + {F_le"}
rhs = TT'*rhs;

% 2nd step: [KK"] = [TT]'[KK_1]
KK    = TT'*KK_1;

KKfree            = KK(free,free);    % stiffness matrix with free indices
GCOORD_temp(free) = KKfree\rhs(free); % solve for the free indices (unknows): 
                                      % {X"(free)} = [KK"(free)]^(-1){F_le"(free)}
GCOORD_temp       = TT*GCOORD_temp; % recover global coordinates {X} = [TT]{X"}

force             = KK_0*GCOORD_temp - F_le; % force vector computed once we know the 
                                             % new positions, {F} = [KK]{X} - {F_le}
GCOORD_temp       = reshape(GCOORD_temp,3,[])';
GCOORD_temp_plot  = GCOORD_temp; % data for plotting

%==========================================================================
% PROJECT BND NODES ONTO BOUNDARIES
%==========================================================================
if r_int > 0 && ~isempty(MESH.pbnd1)
    pbnd1_temp          = GCOORD_temp(nfix+1:nfix+nbnd1,:);    % take bnd1 nodes
    pbnd1_temp_sph      = cartesian2spherical(pbnd1_temp);     % convert to spherical coordinates
    pbnd1_temp_sph(:,3) = r_int*ones(size(pbnd1_temp,1),1);    % change the r-coordinate by r_int (projecting onto bnd 1)
    pbnd1_projected     = spherical2cartesian(pbnd1_temp_sph); % convert to Cartesian coordinates
    GCOORD_temp(nfix+1:nfix+nbnd1,:)             = pbnd1_projected; % substitute pbnd1_projected in GCOORD_temp
end

pbnd2_temp          = GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:); % take bnd2 nodes
pbnd2_temp_sph      = cartesian2spherical(pbnd2_temp);     % convert to spherical coordinates
pbnd2_temp_sph(:,3) = r_ext*ones(size(pbnd2_temp,1),1);    % change the r-coordinate by r_ext (projecting onto bnd 2)
pbnd2_projected     = spherical2cartesian(pbnd2_temp_sph); % convert to Cartesian coordinates
GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:) = pbnd2_projected; % substitute pbnd2_projected in GCOORD_temp

pbnd_ref_face_bot_temp          = GCOORD_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,:); % take pbnd_ref_face_bot
pbnd_ref_face_bot_temp_sph      = cartesian2spherical(pbnd_ref_face_bot_temp);                 % convert to spherical coordinates
pbnd_ref_face_bot_temp_sph(:,3) = (r_ext-d_ref)*ones(size(pbnd_ref_face_bot_temp,1),1);        % change the r-coordinate by r_ext - d_ref (projecting)
pbnd_ref_face_bot_projected     = spherical2cartesian(pbnd_ref_face_bot_temp_sph);             % convert to Cartesian coordinates
GCOORD_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,:) = pbnd_ref_face_bot_projected; % substitute pbnd_ref_face_bot_projected in GCOORD_temp

GCOORD_plot = GCOORD_temp;     % data for plotting

%==========================================================================
% Has been any bnd_ref_face_bot node carried outside of the bottom face of 
% the refined region when we have solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
GCOORD_SPH_temp              = cartesian2spherical(GCOORD_temp);
i_bnd_out_ref_face_bot       = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,1) <= theta_ref_n & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,1) >= theta_ref_s & ...
                               GCOORD_SPH_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,2) <= phi_ref_w   & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,2) >= phi_ref_e; 
n_bnd_out_ref_face_bot       = sum(i_bnd_out_ref_face_bot);
GCOORD_temp_bnd_ref_face_bot = GCOORD_temp(nfix+nbnd1+nbnd2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot,:);
GCOORD_temp_bnd_ref_face_bot(i_bnd_out_ref_face_bot,:) = [];

%==========================================================================
% Has been any bnd_ref_face_1 node carried outside of the face 1 of 
% the refined region when we have solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
i_bnd_out_ref_face_1       = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1,2) <= phi_ref_w     & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1,2) >= phi_ref_e & ...
                             GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1,3) <= r_ext - d_ref & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1,3) >= r_ext; 
n_bnd_out_ref_face_1       = sum(i_bnd_out_ref_face_1);
GCOORD_temp_bnd_ref_face_1 = GCOORD_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1,:);
GCOORD_temp_bnd_ref_face_1(i_bnd_out_ref_face_1,:) = [];

%==========================================================================
% Has been any bnd_ref_face_2 node carried outside of the face 2 of 
% the refined region when we have solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
i_bnd_out_ref_face_2       = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2,1) <= theta_ref_n   & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2,1) >= theta_ref_s & ...
                             GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2,3) <= r_ext - d_ref & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2,3) >= r_ext; 
n_bnd_out_ref_face_2       = sum(i_bnd_out_ref_face_2);
GCOORD_temp_bnd_ref_face_2 = GCOORD_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2,:);
GCOORD_temp_bnd_ref_face_2(i_bnd_out_ref_face_2,:) = [];

%==========================================================================
% Has been any bnd_ref_face_3 node carried outside of the face 3 of 
% the refined region when we have solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
i_bnd_out_ref_face_3       = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3,2) <= phi_ref_w     & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3,2) >= phi_ref_e & ...
                             GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3,3) <= r_ext - d_ref & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3,3) >= r_ext; 
n_bnd_out_ref_face_3       = sum(i_bnd_out_ref_face_3);
GCOORD_temp_bnd_ref_face_3 = GCOORD_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3,:);
GCOORD_temp_bnd_ref_face_3(i_bnd_out_ref_face_3,:) = [];

%==========================================================================
% Has been any bnd_ref_face_4 node carried outside of the face 4 of 
% the refined region when we have solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
i_bnd_out_ref_face_4       = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4,1) <= theta_ref_n   & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4,1) >= theta_ref_s & ...
                             GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4,3) <= r_ext - d_ref & GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4,3) >= r_ext; 
n_bnd_out_ref_face_4       = sum(i_bnd_out_ref_face_4);
GCOORD_temp_bnd_ref_face_4 = GCOORD_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+1:nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4,:);
GCOORD_temp_bnd_ref_face_4(i_bnd_out_ref_face_4,:) = [];

%==========================================================================
% Has been any interior node carried outside of the boundary when we have
% solved for the stiffness matrix? 
% If so, remove them and remove the force in these nodes
%==========================================================================
i_int_out_bnd           = GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4+1:end,3) <= r_int | ...
                          GCOORD_SPH_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4+1:end,3) >= r_ext;
n_int_out_bnd           = sum(i_int_out_bnd);
GCOORD_temp_int_out_bnd = GCOORD_temp(nfix+nbnd1+nbnd2+nbnd_ref_face_bot+nbnd_ref_face_1+nbnd_ref_face_2+nbnd_ref_face_3+nbnd_ref_face_4+1:end,:);
GCOORD_temp_int_out_bnd(i_int_out_bnd,:) = [];

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
if r_int > 0 && ~isempty(MESH.pbnd1)
    MESH.pbnd1  = GCOORD_temp(nfix+1:nfix+nbnd1,:);
end
MESH.pbnd2             = GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:);
MESH.pbnd_ref_face_bot = GCOORD_temp_bnd_ref_face_bot;
MESH.pbnd_ref_face_1   = GCOORD_temp_bnd_ref_face_1;
MESH.pbnd_ref_face_2   = GCOORD_temp_bnd_ref_face_2;
MESH.pbnd_ref_face_3   = GCOORD_temp_bnd_ref_face_3;
MESH.pbnd_ref_face_4   = GCOORD_temp_bnd_ref_face_4;
MESH.pint              = GCOORD_temp_int_out_bnd;
MESH.n_nod_out         = sum([n_bnd_out_ref_face_bot; ...
                              n_bnd_out_ref_face_1;   ...
                              n_bnd_out_ref_face_2;   ...
                              n_bnd_out_ref_face_3;   ...
                              n_bnd_out_ref_face_4;   ...
                              n_int_out_bnd]);
MESH.GCOORD = [MESH.pfix;              ...
               MESH.pbnd1;             ...
               MESH.pbnd2;             ...
               MESH.pbnd_ref_face_bot; ...
               MESH.pbnd_ref_face_1;   ...
               MESH.pbnd_ref_face_2;   ...
               MESH.pbnd_ref_face_3;   ...
               MESH.pbnd_ref_face_4;   ...
               MESH.pint];
MESH.EL2NOD = delaunay(MESH.GCOORD);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_SPH  = cartesian2spherical(MESH.GCOORD);
MESH.EL2NOD = MESH.EL2NOD(~(sum(ismember(MESH.EL2NOD,find(abs(GCOORD_SPH(:,3)-r_int) < 1e-8)),2)==4),:);
%---------------------------------what the line above does--------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_SPH(:,3)-r_int) < 1e-8);
% elements_with_4_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 4;
% remove_elements_with_4_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_4_nodes_on_bnd1),:);
%-----------------------------------------------------------------------------------------------------

force                  = reshape(force,3,[])';
i_fix_bnd1_bnd2_out    = false(nfix+nbnd1+nbnd2,1); % fixed and boundary nodes are not outside the boundaries
force([i_fix_bnd1_bnd2_out;    ...
       i_bnd_out_ref_face_bot; ...
       i_bnd_out_ref_face_1;   ...
       i_bnd_out_ref_face_2;   ...
       i_bnd_out_ref_face_3;   ...
       i_bnd_out_ref_face_4;   ...
       i_int_out_bnd],:)       = []; % remove the force of the outside nodes

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(3)
    subplot(1,2,2)
        simpplot(MESH.GCOORD,MESH.EL2NOD,'p(:,1)<0')
    view(142.5,30)
    text(r_ext/2,-r_ext,r_ext,['iter = ',num2str(MESH.iter)], ...
        'HorizontalAlignment','center','BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
    text(r_ext/2,-r_ext,1.5*r_ext,('Bars (connectivity) after solving nodal positions'), ...
        'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor','black');
    
    if size(MESH.GCOORD,1) < 1500
        figure(35)
        clf
        % initial locations in white
        scatter3(GCOORD_init(:,1),GCOORD_init(:,2),GCOORD_init(:,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1])
        hold on
        if r_int > 0 && ~isempty(MESH.pbnd1)
            % location on the tangent plane to boundary 1 in red (only bnd nodes)
            scatter3(GCOORD_temp_plot(nfix+1:nfix+nbnd1,1), ...
                     GCOORD_temp_plot(nfix+1:nfix+nbnd1,2), ...
                     GCOORD_temp_plot(nfix+1:nfix+nbnd1,3), ...
                     'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
        end
        % location on the tangent plane to boundary 2 in yellow (only bnd nodes)
        scatter3(GCOORD_temp_plot(nfix+nbnd1+1:nfix+nbnd1+nbnd2,1), ...
                 GCOORD_temp_plot(nfix+nbnd1+1:nfix+nbnd1+nbnd2,2), ...
                 GCOORD_temp_plot(nfix+nbnd1+1:nfix+nbnd1+nbnd2,3), ...
                 'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0])
        % final location on surface in blue
        scatter3(GCOORD_plot(:,1),GCOORD_plot(:,2),GCOORD_plot(:,3), ...
                 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
        % nodes outside the boundary in green
        scatter3(GCOORD_plot([temp1;temp2],1),GCOORD_plot([temp1;temp2],2),GCOORD_plot([temp1;temp2],3), ...
                 'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
        line([GCOORD_init(:,1) GCOORD_temp_plot(:,1) GCOORD_plot(:,1)]', ...
             [GCOORD_init(:,2) GCOORD_temp_plot(:,2) GCOORD_plot(:,2)]', ...
             [GCOORD_init(:,3) GCOORD_temp_plot(:,3) GCOORD_plot(:,3)]','Color',[1 0 0])
        title('nodal motion after solve positions')
        text(r_ext,-r_ext,r_ext,['iter = ',num2str(MESH.iter)], ...
            'HorizontalAlignment','center','BackgroundColor',[1 1 0],'Margin',5,'EdgeColor','black');
        view(142.5,30)
        axis equal
        xlabel('X');ylabel('Y');zlabel('Z');
        grid on
        if sum([temp1;temp2]) == 0
            if r_int > 0 && ~isempty(MESH.pbnd1)
                legend('initial location','tangent location boundary 1','tangent location boundary 2','final location')
            else
                legend('initial location','tangent location boundary 2','final location')
            end
        else
            if r_int > 0 && ~isempty(MESH.pbnd1)
                legend('initial location',            ...
                       'tangent location boundary 1', ...
                       'tangent location boundary 2', ...
                       'final location',              ...
                       'nodes outside the boundary')
            else
                legend('initial location',            ...
                       'tangent location boundary 2', ...
                       'final location',              ...
                       'nodes outside the boundary')
            end
        end
        lightGrey = 0.90*[1 1 1];
        [x_sphere,y_sphere,z_sphere] = sphere(50);
        if r_int > 0 && ~isempty(MESH.pbnd1)
            surface(x_sphere*SETTINGS.r_int, ...
                    y_sphere*SETTINGS.r_int, ...
                    z_sphere*SETTINGS.r_int,'FaceColor','none','EdgeColor',lightGrey)
        end
        surface(x_sphere*SETTINGS.r_ext, ...
                y_sphere*SETTINGS.r_ext, ...
                z_sphere*SETTINGS.r_ext,'FaceColor','none','EdgeColor',lightGrey)
    end
end
end % END OF FUNCTION projection_v2