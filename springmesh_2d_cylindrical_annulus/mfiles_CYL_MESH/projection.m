function [MESH] = projection(MESH,STIFFNESS,F_le,SETTINGS)
% Usage: [MESH] = projection(MESH,STIFFNESS,F_le,SETTINGS)
%
% Purpose:
%   Fix boundary nodes motion on the tangent line to the boundary at each
%   node and solve for the answer.
%   We rotate bnd nodes to their own axis domain (X",Y"), while the tangent
%   line is along X". So, by fixing Y" direction displacements, we fix the
%   bnd nodes on their tangent lines to the circumference.
%   We design transformation matrix [TT] for this goal which is in function
%   of theta (colatitude, measured from +Y axis in clockwise direction 
%   (range 0 to 360).
%
%               Y
%               |        Y"
%               |  th   /
%               |___   /
%               |   \ /
%               |    *
%               |   /
%               |  / r
%               | /
%        -------.--------------> X
%               | ---___
%               |       ---___
%               |             -- X"
%               |
%
%   The transformation matrix [TT] to change from local coordinates
%   (x",y") to global coordinates (x,y) is given by:
%   
%   [x]   [  cos(theta)  sin(theta) ][x"]
%   | | = |                         ||  |
%   [y]   [ -sin(theta)  cos(theta) ][y"] (= r for bnd nodes)<-constrain BC
%
% Input:
%   MESH      : [structure]     : structure containing the mesh
%   STIFFNESS : [structure]     : structure containing matrix parameters
%   F_le      : [column vector] : force vector term due to the effect
%                                 of a preferred zero-deformation element
%                                 length
%   SETTINGS  : [structure]     : structure containing mesh settings
%
% Output:
%   MESH      : [structure]     : structure containing the mesh
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
r_int       = SETTINGS.r_int;
r_ext       = SETTINGS.r_ext;
GCOORD      = MESH.GCOORD;
nnod        = size(GCOORD,1);
nfix        = size(MESH.pfix,1);
nbnd1       = size(MESH.pbnd1,1);
nbnd2       = size(MESH.pbnd2,1);
GCOORD_init = GCOORD; % save for plotting

%==========================================================================
% ASSEMBLE TRANSFORMATION MATRIX [TT]
%==========================================================================
% theta angles
if r_int > 0 && ~isempty(MESH.pbnd1)
    theta1 = atan2d(MESH.pbnd1(:,1),MESH.pbnd1(:,2));
end
theta2 = atan2d(MESH.pbnd2(:,1),MESH.pbnd2(:,2));

% pre-allocate memory: number of non-zero slots in [TT]
nnz_T          = 2*nnod + 2*nbnd1 + 2*nbnd2; % = 2*nfix + 4*nbnd1 + 4*nbnd2 + 2*nint
KKi_T          = zeros(nnz_T,1);
KKj_T          = KKi_T;
KK1_T          = ones(nnz_T,1);
KKi_T(1:2*nnod) = (1:2*nnod)';
KKj_T(1:2*nnod) = (1:2*nnod)';

if r_int > 0 && ~isempty(MESH.pbnd1)
    % diagonal slots for boundary 1: cos(theta1) for each node in bnd 1
    % (2 slots for each bnd node).
    KK1_T(2*nfix+1:2:2*(nfix+nbnd1)-1)               = cosd(theta1); % element (1,1) in [TT] for bnd 1 nodes
    KK1_T(2*nfix+2:2:2*(nfix+nbnd1)  )               = cosd(theta1); % element (2,2) in [TT] for bnd 1 nodes
end

% diagonal slots for boundary 2: cos(theta2) for each node in bnd 2
% (2 slots for each bnd node). This will use up 2*nnod slots
KK1_T(2*(nfix+nbnd1)+1:2:2*(nfix+nbnd1+nbnd2)-1) = cosd(theta2); % element (1,1) in [TT] for bnd 2 nodes
KK1_T(2*(nfix+nbnd1)+2:2:2*(nfix+nbnd1+nbnd2)  ) = cosd(theta2); % element (2,2) in [TT] for bnd 2 nodes

if r_int > 0 && ~isempty(MESH.pbnd1)
    % off-diagonal slots for boundary 1. This will use up 2*nbnd1 slots
    KKi_T(2*nnod+1:2*nnod+nbnd1)                     = (2*nfix+1:2:2*(nfix+nbnd1)-1);
    KKj_T(2*nnod+1:2*nnod+nbnd1)                     = (2*nfix+2:2:2*(nfix+nbnd1)  );
    KK1_T(2*nnod+1:2*nnod+nbnd1)                     = sind(theta1);  % element (1,2) in [TT] for bnd 1 nodes
    KKi_T(2*nnod+nbnd1+1:2*nnod+2*nbnd1)             = (2*nfix+2:2:2*(nfix+nbnd1)  );
    KKj_T(2*nnod+nbnd1+1:2*nnod+2*nbnd1)             = (2*nfix+1:2:2*(nfix+nbnd1)-1);
    KK1_T(2*nnod+nbnd1+1:2*nnod+2*nbnd1)             = -sind(theta1); % element (2,1) in [TT] for bnd 1 nodes
end

% off-diagonal slots for boundary 2. This will use up 2*nbnd2 slots
KKi_T(2*nnod+2*nbnd1+1:2*nnod+2*nbnd1+nbnd2)     = (2*(nfix+nbnd1)+1:2:2*(nfix+nbnd1+nbnd2)-1);
KKj_T(2*nnod+2*nbnd1+1:2*nnod+2*nbnd1+nbnd2)     = (2*(nfix+nbnd1)+2:2:2*(nfix+nbnd1+nbnd2)  );
KK1_T(2*nnod+2*nbnd1+1:2*nnod+2*nbnd1+nbnd2)     = sind(theta2);  % element (1,2) in [TT] for bnd 2 nodes
KKi_T(2*nnod+2*nbnd1+nbnd2+1:end)                = (2*(nfix+nbnd1)+2:2:2*(nfix+nbnd1+nbnd2)  );
KKj_T(2*nnod+2*nbnd1+nbnd2+1:end)                = (2*(nfix+nbnd1)+1:2:2*(nfix+nbnd1+nbnd2)-1);
KK1_T(2*nnod+2*nbnd1+nbnd2+1:end)                = -sind(theta2); % element (2,1) in [TT] for bnd 2 nodes

TT = sparse2(KKi_T(:),KKj_T(:),KK1_T(:)); % transformation matrix

%==========================================================================
% ASSEMBLE THE GLOBAL STIFFNESS MATRIX [KK] AND SOLVE FOR ANSWER {X}
%==========================================================================
% The system of equations to be solved before applying the circular BCs is
% given by:
%
% [KK]{X} = {F} + {F_le} where:
%   {F} = 0 since we are solving for steady state of the springs
%   {F_le} is the force vector term due to the effect of a preferred
%    zero-deformation element length
%
% After applying the transformation matrix [TT] (to impose the circular
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

KKi = [STIFFNESS.KKi(:);STIFFNESS.KKi_CrossBars(:)];
KKj = [STIFFNESS.KKj(:);STIFFNESS.KKj_CrossBars(:)];
KK1 = [STIFFNESS.KK1(:);STIFFNESS.KK1_CrossBars(:)];

KK_0 = sparse2(KKi(:),KKj(:),KK1(:)); % MATLAB's sparse matrix maker accumulates all element KK's
free          = [2*nfix+1:2:2*(nfix+nbnd1+nbnd2)-1,...
                 2*(nfix+nbnd1+nbnd2)+1:2*nnod]; % free indices
all_ind       = 1:2*nnod;
all_ind(free) = [];
fix           = all_ind; % fix indices

% Apply transformation matrix (to impose circular BCs) to stiffness
% matrix [KK"] = [TT]'[KK][TT] in two steps:
% 1st step: [KK_1] = [KK_0][TT]
KK_1  = KK_0*TT;
% 2nd step: [KK"] = [TT]'[KK_1]
KK    = TT'*KK_1;
% Apply transformation matrix (transposed) to global coordinates {X"}=[TT]'{X}
% (to impose circular BCs in local coordinates)
GCOORD_temp = GCOORD';
GCOORD_temp = GCOORD_temp(:);
GCOORD_temp = TT'*GCOORD_temp;
if r_int > 0 && ~isempty(MESH.pbnd1)
    % Impose BCs for boundary 1, fixing the y" coordinate (of bnd nodes) with r_int
    GCOORD_temp(2*nfix+2:2:2*(nfix+nbnd1))               = r_int;
end
% Impose BCs for boundary 2, fixing the y" coordinate (of bnd nodes) with r_ext
GCOORD_temp(2*(nfix+nbnd1)+2:2:2*(nfix+nbnd1+nbnd2)) = r_ext;
% Pass imposed BCs to the right hand side in two steps:
% 1st step: rhs = -([KK_0][TT])[TT]'{X} + {F_le}
rhs = -KK_1(:,fix)*GCOORD_temp(fix) + F_le ; % = 0 - KK*pfix
% 2nd step: rhs = [TT]'(-([KK_0][TT])[TT]'{X} + {F_le}) =
%
%               = -([TT]'[KK_0][TT])[TT]'{X} + [TT]'{F_le} =
%
%               = - [KK"(fix)]{X"(fix)} + {F_le"}
rhs = TT'*rhs;

KKfree            = KK(free,free);    % stiffness matrix with free indices
GCOORD_temp(free) = KKfree\rhs(free); % solve for the free indices (unknows): 
                                      % {X"(free)} = [KK"(free)]^(-1){F_le"(free)}
GCOORD_temp       = TT*GCOORD_temp; % recover global coordinates {X} = [TT]{X"}

force       = KK_0*GCOORD_temp - F_le; % force vector computed once we know the 
                                       % new positions, {F} = [KK]{X} - {F_le}
GCOORD_temp = reshape(GCOORD_temp,2,[])';

%==========================================================================
% PROJECT BND NODES ONTO BOUNDARIES
%==========================================================================
GCOORD_temp_plot = GCOORD_temp; % data for plotting

if r_int > 0 && ~isempty(MESH.pbnd1)
    pbnd1 = GCOORD_temp(nfix+1:nfix+nbnd1,:);        % take bnd1 nodes
    theta_temp_bnd1 = atan2d(pbnd1(:,1),pbnd1(:,2)); % theta for new positions of bnd1 nodes
    pbnd1_circ(:,1) = r_int*sind(theta_temp_bnd1);   % x-coordinate after projecting onto the bnd1 (r_int)
    pbnd1_circ(:,2) = r_int*cosd(theta_temp_bnd1);   % y-coordinate after projecting onto the bnd1 (r_int)
    GCOORD_temp(nfix+1:nfix+nbnd1,:)             = pbnd1_circ; % carry points on the tangent line to the boundary 1
end

pbnd2 = GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:); % take bnd2 nodes
theta_temp_bnd2 = atan2d(pbnd2(:,1),pbnd2(:,2)); % theta for new positions of bnd1 nodes
pbnd2_circ(:,1) = r_ext*sind(theta_temp_bnd2);   % x-coordinate after projecting onto the bnd2 (r_ext)
pbnd2_circ(:,2) = r_ext*cosd(theta_temp_bnd2);   % y-coordinate after projecting onto the bnd2 (r_ext)
GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:) = pbnd2_circ; % carry points on the tangent line to the boundary 2

GCOORD_plot = GCOORD_temp;     % data for plotting

%==========================================================================
% Has been any INTERIOR node carried outside of the boundary when we have
% solved for the stiffness matrix? In positive case --> remove them and
% remove the force in these nodes
%==========================================================================
temp1                  = false(nfix+nbnd1+nbnd2,1); % fixed and boundary nodes are not outside the boundaries
r_temp                 = sqrt(GCOORD_temp(nfix+nbnd1+nbnd2+1:end,1).^2 + ...
                              GCOORD_temp(nfix+nbnd1+nbnd2+1:end,2).^2);
temp2                  = r_temp < r_int | r_temp > r_ext;
if strcmp(SETTINGS.mesh,'axisym')
    temp2 = r_temp < r_int | r_temp > r_ext | GCOORD_temp(nfix+nbnd1+nbnd2+1:end,1) < 0;
end
GCOORD_temp([temp1;temp2],:) = [];                  % remove the outside nodes
MESH.n_nod_out         = sum([temp1;temp2]); % number of nodes outside boundaries
force                  = reshape(force,2,[])';
force([temp1;temp2],:) = [];                 % remove the force of the outside nodes

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
if r_int > 0 && ~isempty(MESH.pbnd1)
    MESH.pbnd1  = GCOORD_temp(nfix+1:nfix+nbnd1,:);
end
MESH.pbnd2  = GCOORD_temp(nfix+nbnd1+1:nfix+nbnd1+nbnd2,:);
MESH.pint   = GCOORD_temp(nfix+nbnd1+nbnd2+1:end,:);
MESH.GCOORD = GCOORD_temp;
EL2NOD_temp = delaunay(GCOORD_temp);
% Remove elements created inside the interior boundary (boundary 1)
GCOORD_POL  = cartesian2polar(MESH.GCOORD);
MESH.EL2NOD = EL2NOD_temp(~(sum(ismember(EL2NOD_temp,find(abs(GCOORD_POL(:,2)-r_int) < 1e-8)),2)==3),:);
%---------------------------------what the line above does----------------------------------------------
% nodes_on_bnd1                                  = find(abs(GCOORD_POL(:,2)-r_int) < 1e-8);
% elements_with_3_nodes_on_bnd1                  = sum(ismember(EL2NOD,nodes_on_bnd1),2) == 3;
% remove_elements_with_3_nodes_on_bnd1_in_EL2NOD = EL2NOD(~(elements_with_3_nodes_on_bnd1),:);
%-------------------------------------------------------------------------------------------------------

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(3)
    subplot(1,2,2)
    hold on
    % initial locations in red
    plot(GCOORD_init(:,1),GCOORD_init(:,2),'rx')
    if r_int > 0 && ~isempty(MESH.pbnd1)
        % location on the tangent line to boundary 1 in blue (only bnd nodes)
        plot(GCOORD_temp_plot(nfix+1:nfix+nbnd1,1),...
            GCOORD_temp_plot(nfix+1:nfix+nbnd1,2),'bo')
    end
    % location on the tangent line to boundary 2 in blue (only bnd nodes) 
    plot(GCOORD_temp_plot(nfix+nbnd1+1:nfix+nbnd1+nbnd2,1),...
         GCOORD_temp_plot(nfix+nbnd1+1:nfix+nbnd1+nbnd2,2),'bo')
    % final location on circumference in black
    plot(GCOORD_plot(:,1),GCOORD_plot(:,2),'k+');
    % nodes outside the boundary
    plot(GCOORD_plot([temp1;temp2],1),GCOORD_plot([temp1;temp2],2),'g*')
    line([GCOORD_init(:,1) GCOORD_temp_plot(:,1) GCOORD_plot(:,1)]',...
         [GCOORD_init(:,2) GCOORD_temp_plot(:,2) GCOORD_plot(:,2)]','Color',[1 0 0])
    if sum([temp1;temp2]) == 0
        if r_int > 0 && ~isempty(MESH.pbnd1)
            legend('initial location',...
                   'tangent location boundary 1',...
                   'tangent location boundary 2',...
                   'final location')
        else
            legend('initial location',...
                   'tangent location boundary 2',...
                   'final location')
        end
    else
        if r_int > 0 && ~isempty(MESH.pbnd1)
            legend('initial location',...
                   'tangent location boundary 1',...
                   'tangent location boundary 2',...
                   'final location',...
                   'nodes outside the boundary')
        else
            legend('initial location',...
                   'tangent location boundary 2',...
                   'final location',...
                   'nodes outside the boundary')
        end
    end
    trimesh(MESH.EL2NOD,MESH.GCOORD(:,1),MESH.GCOORD(:,2),'Color',[0 0 0]);
    xlabel('X (km)')
    ylabel('Y (km)')
    axis equal
    axis([-r_ext r_ext -r_ext r_ext])
    title('nodal movement and final bars (connectivity)')
end
end % END OF FUNCTION bc_projection