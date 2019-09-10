function [MESH] = projection(MESH,STIFFNESS,F_le,SETTINGS)
% Usage: [MESH] = projection(MESH,STIFFNESS,F_le,SETTINGS)
%
% Purpose:
%   Fix boundary nodes motion on the boundary segment at each node and
%   solve for the answer.
%   We rotate bnd nodes to their own axis domain (X",Y"), while the segment
%   is along X". So, by fixing Y" direction displacements, we fix the
%   bnd nodes on their segments.
%   We design transformation matrix [TT] for this goal.
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
%   [x]   [  b  a ][x"]
%   | | = |       ||  |
%   [y]   [ -a  b ][y"] (= 0 for bnd nodes)<-constrain BC
%
%   where [a,b] is the slope of the segment (segBC(:,[5 6])
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
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth       = SETTINGS.depth;  % domain depth (km)
length      = SETTINGS.length; % domain length (km)
x0          = SETTINGS.x0;     % point around which the domain is created
z0          = SETTINGS.z0;     % point around which the domain is created
GCOORD      = MESH.GCOORD;
nnod        = size(GCOORD,1);
nfix        = size(MESH.pfix,1);
nbnd        = size(MESH.pbnd,1);
segBC       = MESH.segBC;
pseg        = MESH.pseg;
GCOORD_init = GCOORD; % save for plotting

%==========================================================================
% ASSEMBLE TRANSFORMATION MATRIX [TT]
%==========================================================================

% pre-allocate memory: number of non-zero slots in [TT]
nnz_T           = 2*nnod + 2*nbnd; % = 2*nfix + 4*nbnd + 2*nint
KKi_T           = zeros(nnz_T,1);
KKj_T           = KKi_T;
KK1_T           = ones(nnz_T,1);
KKi_T(1:2*nnod) = (1:2*nnod)';
KKj_T(1:2*nnod) = (1:2*nnod)';

% diagonal slots: b (from segBC) for bnd slots, 1 for all other slots. This will use up the first 2*nnod slots
KK1_T(2*nfix+1:2:2*(nfix+nbnd)-1) = segBC(pseg,6); % element (1,1) in [TT] for bnd nodes
KK1_T(2*nfix+2:2:2*(nfix+nbnd))   = segBC(pseg,6); % element (2,2) in [TT] for bnd nodes

% off-diagonal slots: a and -a (from segBC). This will use up the last 2*nnod slots 
KKi_T(2*nnod+1:2*nnod+nbnd)       = (2*nfix+1:2:2*(nfix+nbnd)-1);
KKj_T(2*nnod+1:2*nnod+nbnd)       = (2*nfix+2:2:2*(nfix+nbnd)  );
KK1_T(2*nnod+1:2*nnod+nbnd)       = segBC(pseg,5); % element (1,2) in [TT] for bnd nodes

KKi_T(2*nnod+nbnd+1:end)          = (2*nfix+2:2:2*(nfix+nbnd)  );
KKj_T(2*nnod+nbnd+1:end)          = (2*nfix+1:2:2*(nfix+nbnd)-1);
KK1_T(2*nnod+nbnd+1:end)          = -segBC(pseg,5); % element (2,1) in [TT] for bnd nodes

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
free          = [2*nfix+1:2:2*(nfix+nbnd)-1,...
                 2*(nfix+nbnd)+1:2*nnod]; % free indices
all_ind       = 1:2*nnod;
all_ind(free) = [];
fix           = all_ind; % fix indices

% Apply transformation matrix (to impose BCs) to stiffness matrix [KK"] = [TT]'[KK][TT] in two steps:
% 1st step: [KK_1] = [KK_0][TT]
KK_1  = KK_0*TT;
% 2nd step: [KK"] = [TT]'[KK_1]
KK    = TT'*KK_1;
% Apply transformation matrix (transposed) to global coordinates {X"}=[TT]'{X} (to impose BCs in local coordinates)
GCOORD_temp = GCOORD';
GCOORD_temp = GCOORD_temp(:);
GCOORD_temp = TT'*GCOORD_temp;
% Impose BCs for boundaries, fixing the y" coordinate (of bnd nodes)
GCOORD_temp(2*nfix+2:2:2*(nfix+nbnd)) = (segBC(pseg,5).*(segBC(pseg,1)) + segBC(pseg,6).*(segBC(pseg,2)));

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

GCOORD_plot = GCOORD_temp;     % data for plotting

%==========================================================================
% Has been any INTERIOR node carried outside of the boundary when we have
% solved for the stiffness matrix? In positive case --> remove them and
% remove the force in these nodes
%==========================================================================
xmin                   = x0 - length/2;
xmax                   = x0 + length/2;
zmin                   = z0 - depth;
zmax                   = z0;
temp1                  = false(nfix+nbnd,1); % fixed and boundary nodes are not outside the boundaries
x_temp                 = GCOORD_temp(nfix+nbnd+1:end,1);
z_temp                 = GCOORD_temp(nfix+nbnd+1:end,2);
temp2                  = x_temp > xmax | x_temp < xmin | z_temp > zmax | z_temp < zmin;
GCOORD_temp([temp1;temp2],:) = [];                  % remove the outside nodes
MESH.n_nod_out         = sum([temp1;temp2]); % number of nodes outside boundaries
force                  = reshape(force,2,[])';
force([temp1;temp2],:) = [];                 % remove the force of the outside nodes

%==========================================================================
% DATA FOR OUTPUT
%==========================================================================
MESH.pbnd   = GCOORD_temp(nfix+1:nfix+nbnd,:);
MESH.pint   = GCOORD_temp(nfix+nbnd+1:end,:);
MESH.GCOORD = GCOORD_temp;
MESH.EL2NOD = delaunay(GCOORD_temp);

%==========================================================================
% PLOTS
%==========================================================================
if SETTINGS.show_figs
    figure(3)
    subplot(2,1,2)
    hold on
    % initial locations in red
    plot(GCOORD_init(:,1),GCOORD_init(:,2),'rx')
    % final location in black
    plot(GCOORD_plot(:,1),GCOORD_plot(:,2),'k+');
    % nodes outside the boundary
    plot(GCOORD_plot([temp1;temp2],1),GCOORD_plot([temp1;temp2],2),'g*')
    line([GCOORD_init(:,1) GCOORD_plot(:,1)]',...
         [GCOORD_init(:,2) GCOORD_plot(:,2)]','Color',[1 0 0])
    legend('initial location',...
           'final location',...
           'nodes outside the boundary')
    trimesh(MESH.EL2NOD,MESH.GCOORD(:,1),MESH.GCOORD(:,2),'Color',[0 0 0]);
    xlabel('X (km)')
    ylabel('Y (km)')
    axis equal
    title('nodal movement and final bars (connectivity)')
end
end % END OF FUNCTION bc_projection