function [MESH,STIFFNESS,F_le] = assembly_stiffness(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: [MESH,STIFFNESS,F_le] = assembly_stiffness(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   Assembly the stifnness matrix and the force vector.
%   Stifness matrix (6x6) for all bars is given by:
%
%          [ ca*cb   0   ]                                                          
%          |             |                                                          
%          | ca*sb   0   |                                                          
%          |             |                                                          
%          |  sa     0   |[ -1    1 ][ ca*cb  ca*sb   sa      0      0      0 ]
%   [KK] = |             ||         ||                                        |
%          |   0   ca*cb |[  1   -1 ][   0      0      0    ca*cb  ca*sb   sa ]
%          |             |                                                          
%          |   0   ca*sb |                                                          
%          |             |                                                          
%          [   0    sa   ]                                                          
%   
%          [ -ca2cb2    -ca2sbcb    -sacacb    ca2cb2    ca2sbcb    sacacb ]
%          |                                                               |
%          | -ca2sbcb   -ca2sb2     -sacasb    ca2sbcb   ca2sb2     sacasb |
%          |                                                               |
%          | -sacacb    -sacasb     -sa2       sacacb    sacasb     sa2    |
%        = |                                                               |
%          |  ca2cb2     ca2sbcb     sacacb   -ca2cb2   -ca2sbcb   -sacacb |
%          |                                                               |
%          |  ca2sbcb    ca2sb2      sacasb   -ca2sbcb  -ca2sb2    -sacasb |
%          |                                                               |
%          [  sacacb     sacasb      sa2      -sacacb   -sacasb    -sa2    ]       
%   
%   Force vector term (6x1) due to the effect of a preferred 
%   zero-deformation element length for all bars is given by:
%   
%            [ ca*cb   0   ]                       [  ca*cb  ]
%            |             |                       |         |
%            | ca*sb   0   |                       |  ca*sb  |
%            |             |                       |         |
%            |  sa     0   |[ -1    1 ][ 0  ]      |   sa    |
%   {F_le} = |             ||         ||    | = L0*|         |
%            |   0   ca*cb |[  1   -1 ][ L0 ]      | -ca*cb  |
%            |             |                       |         |
%            |   0   ca*sb |                       | -ca*sb  |
%            |             |                       |         |
%            [   0    sa   ]                       [  -sa    ]
%   
%   where sa,ca,sb,cb are sin(alpha),cos(alpha),sin(beta),cos(beta)
%   respectively. alpha and beta are analogous to latitude and longitude
%   respectively.
%
% Input:
%   MESH       : [structure]     : structure containing the mesh
%   GUIDE_MESH : [structure]     : structure containing guide mesh
%   INTERFACE  : [structure]     : structure containing interface settings
%   SETTINGS   : [structure]     : structure containing mesh settings
%
% Output:
%   MESH       : [structure]     : structure containing the mesh
%   STIFFNESS  : [structure]     : structure containing matrix parameters
%   F_le       : [column vector] : force vector term due to the effect of 
%                                  a preferred zero-deformation element
%                                  length
%
% JMT Jun 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% GEOMETRY OF BARS/SPRINGS (DIRECTION, NEUTRAL LENGTH)
%==========================================================================
[MESH.L0,~,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS);
alpha   = atan2d(barvec(:,3),sqrt(barvec(:,1).^2+barvec(:,2).^2)); % analogous to latitude (°)
beta    = atan2d(barvec(:,2),barvec(:,1));                         % analogous to longitude (°)
sa      = sind(alpha);
ca      = cosd(alpha);
sb      = sind(beta);
cb      = cosd(beta);

% Rotation matrix (2x6) for all bars is given by
%
%        [ ca*cb    ca*sb     sa        0        0        0 ]
% [RR] = |                                                  |
%        [   0        0        0      ca*cb    ca*sb     sa ]
ca2cb2  = (ca.*cb).^2;   % cos(alpha)^2*cos(betha)^2
ca2sb2  = (ca.*sb).^2;   % cos(alpha)^2*sin(betha)^2
sa2     = sa.^2;         % sin(alpha)^2
ca2sbcb = ca.^2.*sb.*cb; % cos(alpha)^2*sin(betha)*cos(betha)
sacacb  = sa.*ca.*cb;    % sin(alpha)*cos(alpha)*cos(betha)
sacasb  = sa.*ca.*sb;    % sin(alpha)*cos(alpha)*sin(betha)

%==========================================================================
% BLOCKING PARAMETERS
%==========================================================================
% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
nelblk        = min(900, nbar);    % how many bars are in one block
nblk          = ceil(nbar/nelblk); % number of blocks
il            = 1;
iu            = nelblk;
KK_this_block = zeros(6,6,nelblk);
f_this_block  = zeros(6,1,nelblk);
STIFFNESS.KK1 = zeros(36*nbar,1);
STIFFNESS.KKi = zeros(36*nbar,1);
STIFFNESS.KKj = zeros(36*nbar,1);
dofs          = zeros(1,6*nelblk);
rows          = zeros(6,1,nelblk);
cols          = zeros(1,6,nelblk); 
F1            = zeros(6*nbar,1);
Fi            = zeros(6*nbar,1);   % preallocate for speed

%==========================================================================
% BLOCK LOOP - STIFFNESS MATRIX COMPUTATION
%==========================================================================
for ib = 1:nblk
    sa_this_block        = sa(il:iu);
    ca_this_block        = ca(il:iu);
    sb_this_block        = sb(il:iu);
    cb_this_block        = cb(il:iu);
    ca2cb2_this_block    = ca2cb2(il:iu);
    ca2sb2_this_block    = ca2sb2(il:iu);
    sa2_this_block       = sa2(il:iu);
    ca2sbcb_this_block   = ca2sbcb(il:iu);
    sacacb_this_block    = sacacb(il:iu);
    sacasb_this_block    = sacasb(il:iu);
    el_this_blo          = il:iu;
    L_0                  = MESH.L0(el_this_blo);  % L_0 is a nelblo-by-1 vector
    
    dofs(1:3:nelblk*6-2) = (3*bars(el_this_blo,:)-2)'; % x-global dofs for this bar element
    dofs(2:3:nelblk*6-1) = (3*bars(el_this_blo,:)-1)'; % y-global dofs for this bar element
    dofs(3:3:nelblk*6)   = 3*bars(el_this_blo,:)';     % z-global dofs for this bar element
    rows(1:3:nelblk*6-2) = (3*bars(el_this_blo,:)-2)'; % x-global dofs for this bar element
    rows(2:3:nelblk*6-1) = (3*bars(el_this_blo,:)-1)'; % y-global dofs for this bar element
    rows(3:3:nelblk*6)   = 3*bars(el_this_blo,:)';     % z-global dofs for this bar element
    cols(1:3:nelblk*6-2) = (3*bars(el_this_blo,:)-2)'; % x-global dofs for this bar element
    cols(2:3:nelblk*6-1) = (3*bars(el_this_blo,:)-1)'; % y-global dofs for this bar element
    cols(3:3:nelblk*6)   = 3*bars(el_this_blo,:)';     % z-global dofs for this bar element
    rows_blo             = repmat(rows,[1,6,1]);
    cols_blo             = repmat(cols,[6,1,1]);
    
    % KK_this_block = RR_this_block' * [-1, 1; 1, -1] * RR_this_block;
    KK_this_block(1,1,:) = -ca2cb2_this_block;
    KK_this_block(1,2,:) = -ca2sbcb_this_block;
    KK_this_block(1,3,:) = -sacacb_this_block;
    KK_this_block(1,4,:) =  ca2cb2_this_block;
    KK_this_block(1,5,:) =  ca2sbcb_this_block;
    KK_this_block(1,6,:) =  sacacb_this_block;
    KK_this_block(2,1,:) = -ca2sbcb_this_block;
    KK_this_block(2,2,:) = -ca2sb2_this_block;
    KK_this_block(2,3,:) = -sacasb_this_block;
    KK_this_block(2,4,:) =  ca2sbcb_this_block;
    KK_this_block(2,5,:) =  ca2sb2_this_block;
    KK_this_block(2,6,:) =  sacasb_this_block;
    KK_this_block(3,1,:) = -sacacb_this_block;
    KK_this_block(3,2,:) = -sacasb_this_block;
    KK_this_block(3,3,:) = -sa2_this_block;
    KK_this_block(3,4,:) =  sacacb_this_block;
    KK_this_block(3,5,:) =  sacasb_this_block;
    KK_this_block(3,6,:) =  sa2_this_block;
    KK_this_block(4,1,:) =  ca2cb2_this_block;
    KK_this_block(4,2,:) =  ca2sbcb_this_block;
    KK_this_block(4,3,:) =  sacacb_this_block;
    KK_this_block(4,4,:) = -ca2cb2_this_block;
    KK_this_block(4,5,:) = -ca2sbcb_this_block;
    KK_this_block(4,6,:) = -sacacb_this_block;
    KK_this_block(5,1,:) =  ca2sbcb_this_block;
    KK_this_block(5,2,:) =  ca2sb2_this_block;
    KK_this_block(5,3,:) =  sacasb_this_block;
    KK_this_block(5,4,:) = -ca2sbcb_this_block;
    KK_this_block(5,5,:) = -ca2sb2_this_block;
    KK_this_block(5,6,:) = -sacasb_this_block;
    KK_this_block(6,1,:) =  sacacb_this_block;
    KK_this_block(6,2,:) =  sacasb_this_block;
    KK_this_block(6,3,:) =  sa2_this_block;
    KK_this_block(6,4,:) = -sacacb_this_block;
    KK_this_block(6,5,:) = -sacasb_this_block;
    KK_this_block(6,6,:) = -sa2_this_block;
    
    STIFFNESS.KK1(36*il-35:36*iu)  = KK_this_block(:);
    STIFFNESS.KKi(36*il-35:36*iu)  = rows_blo(:);
    STIFFNESS.KKj(36*il-35:36*iu)  = cols_blo(:);
    
    % f_this_block = RR_this_block' * [-1, 1; 1, -1] * [0 L0]'
    cacbL_0              = ca_this_block.*cb_this_block.*L_0;
    casbL_0              = ca_this_block.*sb_this_block.*L_0;
    saL_0                = sa_this_block.*L_0;
    f_this_block(1,1,:)  =  cacbL_0;
    f_this_block(2,1,:)  =  casbL_0;
    f_this_block(3,1,:)  =  saL_0;
    f_this_block(4,1,:)  = -cacbL_0;
    f_this_block(5,1,:)  = -casbL_0;
    f_this_block(6,1,:)  = -saL_0;
    
    F1(6*il-5:6*iu)      = f_this_block(:);
    Fi(6*il-5:6*iu)      = dofs;
    il                   = il+nelblk;
    
    % adjust the block size for the last block, because usually it's smaller
    if(ib == nblk-1)
        nelblk        = nbar - iu;
        KK_this_block = zeros(6,6,nelblk);
        f_this_block  = zeros(6,1,nelblk);
        dofs          = zeros(1,6*nelblk);
        rows          = zeros(6,1,nelblk);
        cols          = zeros(1,6,nelblk);
    end
    
    iu = iu+nelblk;
    
end

F_le = accumarray(Fi,F1); % force vector term due to the effect of a 
                          % preferred zero-deformation element length

if SETTINGS.show_figs
    figure(3)
    clf
    subplot(1,2,1)
    simpplot(MESH.GCOORD,MESH.EL2NOD,'p(:,1)<0')
    view(142.5,30)
    text(SETTINGS.r_ext/2,-SETTINGS.r_ext,SETTINGS.r_ext,['iter = ',num2str(MESH.iter)], ...
        'HorizontalAlignment','center','BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
    text(SETTINGS.r_ext/2,-SETTINGS.r_ext,1.5*SETTINGS.r_ext,  ...
         'Bars (connectivity) before solving nodal positions', ...
         'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor','black');
end
end % END OF FUNCTION assembly_stiffness