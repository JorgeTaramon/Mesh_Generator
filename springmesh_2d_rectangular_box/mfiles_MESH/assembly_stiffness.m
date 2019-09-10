function [MESH,STIFFNESS,F_le] = assembly_stiffness(MESH,GUIDE_MESH,SETTINGS)
% Usage: [MESH,STIFFNESS,F_le] = assembly_stiffness(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Assembly the stifnness matrix and the force vector.
%   Stifness matrix (4x4) for all bars is given by:
%   
%          [ c   0 ]                               [ -cc   -sc    cc    sc ]
%          |       |                               |                       |
%          | s   0 |[ -1    1 ][ c   s   0   0 ]   | -sc   -ss    sc    ss | 
%   [KK] = |       ||         ||               | = |                       |
%          | 0   c |[  1   -1 ][ 0   0   c   s ]   |  cc    sc   -cc   -sc |
%          |       |                               |                       |
%          [ 0   s ]                               [  sc    ss   -sc   -ss ]
%
%   Force vector term (4x1) due to the effect of a preferred 
%   zero-deformation element length for all bars is given by:
%   
%            [ c   0 ]                       [  c ]
%            |       |                       |    |
%            | s   0 |[ -1    1 ][ 0  ]      |  s |
%   {F_le} = |       ||         ||    | = L0*|    |
%            | 0   c |[  1   -1 ][ L0 ]      | -c |
%            |       |                       |    |
%            [ 0   s ]                       [ -s ]
%
% Input:
%   MESH       : [structure]     : structure containing the mesh
%   GUIDE_MESH : [structure]     : structure containing guide mesh
%   SETTINGS   : [structure]     : structure containing mesh settings
%
% Output:
%   MESH       : [structure]     : structure containing the mesh
%   STIFFNESS  : [structure]     : structure containing matrix parameters
%   F_le       : [column vector] : force vector term due to the effect of 
%                                  a preferred zero-deformation element
%                                  length
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% GEOMETRY OF BARS/SPRINGS (DIRECTION, NEUTRAL LENGTH)
%==========================================================================
[MESH.L0,L,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,SETTINGS);

% Rotation matrix (2x4) for all bars is given by
%
%        [ c   s   0   0 ]
% [RR] = |               |
%        [ 0   0   c   s ]
s  = barvec(:,2) ./ L(:); % sin(alpha) for RR matrix
c  = barvec(:,1) ./ L(:); % cos(alpha) for RR matrix
cc = c.^2;                % (cos(alpha))^2 for each bar
ss = s.^2;                % (sin(alpha))^2 for each bar
sc = c.*s;                % sin(alpha)*cos(alpha) for each bar

%==========================================================================
% BLOCKING PARAMETERS
%==========================================================================
% MILAMIN method (blocks of elements handled at once)
% see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
nelblk        = min(900, nbar);    % how many bars are in one block
nblk          = ceil(nbar/nelblk); % number of blocks
il            = 1;
iu            = nelblk;
KK_this_block = zeros(4,4,nelblk);
f_this_block  = zeros(4,1,nelblk);
STIFFNESS.KK1 = zeros(16*nbar,1);
STIFFNESS.KKi = zeros(16*nbar,1);
STIFFNESS.KKj = zeros(16*nbar,1);
dofs          = zeros(1,4*nelblk);
rows          = zeros(4,1,nelblk);
cols          = zeros(1,4,nelblk); 
F1            = zeros(4*nbar,1);
Fi            = zeros(4*nbar,1);   % preallocate for speed

%==========================================================================
% BLOCK LOOP - STIFFNESS MATRIX COMPUTATION
%==========================================================================
for ib = 1:nblk
    c_this_block          = c(il:iu);
    s_this_block          = s(il:iu);
    cc_this_block         = cc(il:iu);
    ss_this_block         = ss(il:iu);
    sc_this_block         = sc(il:iu);
    el_this_block         = il:iu;
    L_0                   = MESH.L0(el_this_block); % L_0 is a nelblk-by-1 vector
    
    dofs (1:2:nelblk*4-1) = (2*bars(el_this_block,:)-1)'; % x-global dofs for this bar element
    dofs (2:2:nelblk*4)   =  2*bars(el_this_block,:)';    % y-global dofs for this bar element
    rows (1:2:nelblk*4-1) = (2*bars(el_this_block,:)-1)'; % x-global dofs for this bar element
    rows (2:2:nelblk*4)   =  2*bars(el_this_block,:)';    % y-global dofs for this bar element
    cols (1:2:nelblk*4-1) = (2*bars(el_this_block,:)-1)'; % x-global dofs for this bar element
    cols (2:2:nelblk*4)   =  2*bars(el_this_block,:)';    % y-global dofs for this bar element
    rows_blo              = repmat(rows,[1,4,1]);
    cols_blo              = repmat(cols,[4,1,1]);
    
    % KK_this_block = RR_this_block' * [-1, 1; 1, -1] * RR_this_block;
    KK_this_block(1,1,:)  = -cc_this_block;
    KK_this_block(1,2,:)  = -sc_this_block;
    KK_this_block(1,3,:)  =  cc_this_block;
    KK_this_block(1,4,:)  =  sc_this_block;
    KK_this_block(2,1,:)  = -sc_this_block;
    KK_this_block(2,2,:)  = -ss_this_block;
    KK_this_block(2,3,:)  =  sc_this_block;
    KK_this_block(2,4,:)  =  ss_this_block;
    KK_this_block(3,1,:)  =  cc_this_block;
    KK_this_block(3,2,:)  =  sc_this_block;
    KK_this_block(3,3,:)  = -cc_this_block;
    KK_this_block(3,4,:)  = -sc_this_block;
    KK_this_block(4,1,:)  =  sc_this_block;
    KK_this_block(4,2,:)  =  ss_this_block;
    KK_this_block(4,3,:)  = -sc_this_block;
    KK_this_block(4,4,:)  = -ss_this_block;
    
    STIFFNESS.KK1(16*il-15:16*iu) = KK_this_block(:);
    STIFFNESS.KKi(16*il-15:16*iu) = rows_blo(:);
    STIFFNESS.KKj(16*il-15:16*iu) = cols_blo(:);
    
    % f_this_block = RR_this_block' * [-1, 1; 1, -1] * [0 L0]'
    cL_0                  = c_this_block.*L_0;
    sL_0                  = s_this_block.*L_0;
    f_this_block(1,1,:)   =  cL_0;
    f_this_block(2,1,:)   =  sL_0;
    f_this_block(3,1,:)   = -cL_0;
    f_this_block(4,1,:)   = -sL_0;
    
    F1(4*il-3:4*iu)       = f_this_block(:);
    Fi(4*il-3:4*iu)       = dofs;
    il                    = il+nelblk;
    
    % adjust the block size for the last block, because usually it is smaller
    if(ib == nblk-1)
        nelblk        = nbar - iu;
        KK_this_block = zeros(4,4,nelblk);
        f_this_block  = zeros(4,1,nelblk);
        dofs          = zeros(1,4*nelblk);
        rows          = zeros(4,1,nelblk);
        cols          = zeros(1,4,nelblk);
    end
    
    iu = iu+nelblk;
    
end
F_le = accumarray(Fi,F1); % force vector term due to the effect of a 
                          % preferred zero-deformation element length

if SETTINGS.show_figs
    figure(3)
    clf
    subplot(2,1,1)
    trimesh(MESH.EL2NOD,MESH.GCOORD(:,1),MESH.GCOORD(:,2),'Color',[0 0 0])
    hold on
    plot(MESH.GCOORD(:,1),MESH.GCOORD(:,2),'r.')
    xmin = SETTINGS.x0 - SETTINGS.length/2;
    text(xmin*0.8,-200,['iter = ',num2str(MESH.iter)],'HorizontalAlignment', ...
        'center','BackgroundColor',[1 1 0.5],'Margin',5,'EdgeColor','black');
    xlabel('X (km)')
    ylabel('Y (km)')
    title('Bars (connectivity) before solving nodal positions')
    axis equal
end
end % END OF FUNCTION assembly_stiffness