function [KK1_CrossBars,KKi_CrossBars,KKj_CrossBars,F_le] = cross_bars(MESH,F_le,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: [KK1_CrossBars,KKi_CrossBars,KKj_CrossBars,F_le] = cross_bars(MESH,F_le,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%  We add 4 new bars(springs) in every tetrahedron (i.e., spring 1-5, 2-6,
%  3-7, 4-8), asking them to be the same length, therefore these 4 bars 
%  will help generating a better tetrahedron (closer to a regular one).
%
%                4                          We are going to generate 4 new
%                .                          nodes in the barycentre of each
%               / \\                        side (5,6,7 and 8) but we are 
%              /  .\  \                     not going to add new dofs. The
%             /     \    \                  (x,y,z) for nodes 5,6,7 and 8
%            /     . \      \               are only used to calculate the
%           /         \        \            direction cosines (in function
%          /        .6 \     .5   \         of alpha and beta) for the
%         /             \.           \      rotation matrix. After that all
%        /      7.   . . \.   . __.__ . 3   the information related to 
%       /        .       _\___---   /       nodes 5,6,7 and 8 will be 
%      /     .  _____-.-- .\      /         divided and then added to nodes
%     / _.___---       .8   \   /           1,2,3 and 4. For example: for 
%    .-----------------------./             spring 1-5, the information 
%   1                         2             related to node 1 will be fully
%                                           added to node 1, but the 
%   information related to node 5 will be splited into 3 parts and be added
%   to nodes 2,3 and 4. 
%
% Input:
%   MESH          : [structure]     : structure containing the mesh
%   F_le          : [column vector] : force vector term due to the effect
%                                     of a preferred zero-deformation
%                                     element length
%   GUIDE_MESH    : [structure]     : structure containing guide mesh
%   INTERFACE     : [structure]     : structure containing interface
%                                     settings 
%   SETTINGS      : [structure]     : structure containing mesh settings
%
% Output:
%   KK1_CrossBars : [matrix]        : numerical value for the element (i,j)
%                                     in the stiffness matrix
%   KKi_CrossBars : [matrix]        : i index (arrow) in the stiffness
%                                     matrix 
%   KKj_CrossBars : [matrix]        : j index (column) in the stiffness
%                                     matrix
%   F_le          : [column vector] : force vector term due to the effect 
%                                     of a preferred zero-deformation 
%                                     element length after applying cross
%                                     bars 
% JMT Jan 2016
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

GCOORD   = MESH.GCOORD;
EL2NOD   = MESH.EL2NOD;
k        = 1; % define spring constant (by default is 1)
nel_tet  = length(EL2NOD);
bars     = [EL2NOD(:,[1,2]); ...
            EL2NOD(:,[1,3]); ...
            EL2NOD(:,[1,4]); ...
            EL2NOD(:,[2,3]); ...
            EL2NOD(:,[2,4]); ...
            EL2NOD(:,[3,4])]; % bars defining each tetrahedron
barmid   = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoint of each bar
% compute L0 for all bars (even those bars whose ends are fixed nodes)
switch SETTINGS.refinement
    case 'regular'
        L0 = SETTINGS.h0*ones(size(barmid,1),1);
    case 'interface'
        L0 = bar_L0_interface(barmid,SETTINGS.h0,INTERFACE);
    case 'guide_mesh'
        L0 = bar_L0_guide(barmid,GUIDE_MESH);
end
[~,~,J]  = unique(sort(bars,2),'rows');
% List of bar vector for each cross bar of the triangle
barvec15 = (GCOORD(EL2NOD(:,2),:) + GCOORD(EL2NOD(:,3),:) + GCOORD(EL2NOD(:,4),:))/3 - GCOORD(EL2NOD(:,1),:); % vector from node 1 to node 5 (barycenter of side formed by nodes 2,3,4)
barvec26 = (GCOORD(EL2NOD(:,1),:) + GCOORD(EL2NOD(:,3),:) + GCOORD(EL2NOD(:,4),:))/3 - GCOORD(EL2NOD(:,2),:); % vector from node 2 to node 6 (barycenter of side formed by nodes 1,3,4)
barvec37 = (GCOORD(EL2NOD(:,1),:) + GCOORD(EL2NOD(:,2),:) + GCOORD(EL2NOD(:,4),:))/3 - GCOORD(EL2NOD(:,3),:); % vector from node 3 to node 7 (barycenter of side formed by nodes 1,2,4)
barvec48 = (GCOORD(EL2NOD(:,1),:) + GCOORD(EL2NOD(:,2),:) + GCOORD(EL2NOD(:,3),:))/3 - GCOORD(EL2NOD(:,4),:); % vector from node 4 to node 8 (barycenter of side formed by nodes 1,2,3)
% Generate all the sin(alpha),cos(alpha),sin(betha),cos(betha) for each bar
betha15  = atan2d(barvec15(:,2),barvec15(:,1));
alpha15  = atan2d(barvec15(:,3),sqrt(barvec15(:,1).^2+barvec15(:,2).^2));
sa15     = sind(alpha15);
ca15     = cosd(alpha15);
sb15     = sind(betha15);
cb15     = cosd(betha15);

betha26  = atan2d(barvec26(:,2),barvec26(:,1));
alpha26  = atan2d(barvec26(:,3),sqrt(barvec26(:,1).^2+barvec26(:,2).^2));
sa26     = sind(alpha26);
ca26     = cosd(alpha26);
sb26     = sind(betha26);
cb26     = cosd(betha26);

betha37  = atan2d(barvec37(:,2),barvec37(:,1));
alpha37  = atan2d(barvec37(:,3),sqrt(barvec37(:,1).^2+barvec37(:,2).^2));
sa37     = sind(alpha37);
ca37     = cosd(alpha37);
sb37     = sind(betha37);
cb37     = cosd(betha37);

betha48  = atan2d(barvec48(:,2),barvec48(:,1));
alpha48  = atan2d(barvec48(:,3),sqrt(barvec48(:,1).^2+barvec48(:,2).^2));
sa48     = sind(alpha48);
ca48     = cosd(alpha48);
sb48     = sind(betha48);
cb48     = cosd(betha48);

% Generate all the L0 info: the 4 cross bars in each tetrahedron share the  
% same L0, and this L0 = sqrt(6)/3*average(L0 value for the edges)
L0_CrossBars = sqrt(6)/3*mean(reshape(L0(J),nel_tet,6),2); % nel_tet-by-1
%---------------------------- what the above line does -------------------------------
% L0_with_duplication  = MESH.L0(J);
% L0_reshape           = reshape(L0_with_duplication,nel_tet,6); % nel_tet-by-6 matrix
% L0_mean_for_each_tet = mean(L0_reshape,2);                     % nel_tet-by-1 vector
% L0_CrossBars         = sqrt(6)/3*L0_mean_for_each_tet;
%-------------------------------------------------------------------------------------

KK1_CrossBars  = zeros(144,nel_tet);
KKi_CrossBars  = zeros(144,nel_tet);
KKj_CrossBars  = zeros(144,nel_tet);
dofs_CrossBars = zeros(1,12);

TT15 = [1  0  0  0  0  0  0  0  0  0  0  0  ;...
        0  1  0  0  0  0  0  0  0  0  0  0  ;...
        0  0  1  0  0  0  0  0  0  0  0  0  ;...
        0  0  0 1/3 0  0 1/3 0  0 1/3 0  0  ;...
        0  0  0  0 1/3 0  0 1/3 0  0 1/3 0  ;...
        0  0  0  0  0 1/3 0  0 1/3 0  0 1/3];

TT26 = [0  0  0  1  0  0  0  0  0  0  0  0  ;...
        0  0  0  0  1  0  0  0  0  0  0  0  ;...
        0  0  0  0  0  1  0  0  0  0  0  0  ;...
       1/3 0  0  0  0  0 1/3 0  0 1/3 0  0  ;...
        0 1/3 0  0  0  0  0 1/3 0  0 1/3 0  ;...
        0  0 1/3 0  0  0  0  0 1/3 0  0 1/3];

TT37 = [0  0  0  0  0  0  1  0  0  0  0  0  ;...
        0  0  0  0  0  0  0  1  0  0  0  0  ;...
        0  0  0  0  0  0  0  0  1  0  0  0  ;...
       1/3 0  0 1/3 0  0  0  0  0 1/3 0  0  ;...
        0 1/3 0  0 1/3 0  0  0  0  0 1/3 0  ;...
        0  0 1/3 0  0 1/3 0  0  0  0  0 1/3];

TT48 = [0  0  0  0  0  0  0  0  0  1  0  0 ;...
        0  0  0  0  0  0  0  0  0  0  1  0 ;...
        0  0  0  0  0  0  0  0  0  0  0  1 ;...
       1/3 0  0 1/3 0  0 1/3 0  0  0  0  0 ;...
        0 1/3 0  0 1/3 0  0 1/3 0  0  0  0 ;...
        0  0 1/3 0  0 1/3 0  0 1/3 0  0  0];

% loop over all the tetrahedrons, calculate the fel and KK info, also 
% generates the indices for all the fel, KK values (for sparse matrix)
for iel = 1:nel_tet
    
    ca2cb215  = ca15(iel)^2*cb15(iel)^2;
    ca2sb215  = ca15(iel)^2*sb15(iel)^2;
    sa215     = sa15(iel)^2;
    ca2sbcb15 = ca15(iel)^2*sb15(iel)*cb15(iel);
    sacacb15  = sa15(iel)*ca15(iel)*cb15(iel);
    sacasb15  = sa15(iel)*ca15(iel)*sb15(iel);
    
    KK15      = k*[-ca2cb215  -ca2sbcb15 -sacacb15  ca2cb215   ca2sbcb15  sacacb15 ;...
                   -ca2sbcb15 -ca2sb215  -sacasb15  ca2sbcb15  ca2sb215   sacasb15 ;...
                   -sacacb15  -sacasb15  -sa215     sacacb15   sacasb15   sa215    ;...
                    ca2cb215   ca2sbcb15  sacacb15 -ca2cb215  -ca2sbcb15 -sacacb15 ;...
                    ca2sbcb15  ca2sb215   sacasb15 -ca2sbcb15 -ca2sb215  -sacasb15 ;...
                    sacacb15   sacasb15   sa215    -sacacb15  -sacasb15  -sa215   ];  % Stiffness matrix for bar vector 1-5
    
    ca2cb226  = ca26(iel)^2*cb26(iel)^2;
    ca2sb226  = ca26(iel)^2*sb26(iel)^2;
    sa226     = sa26(iel)^2;
    ca2sbcb26 = ca26(iel)^2*sb26(iel)*cb26(iel);
    sacacb26  = sa26(iel)*ca26(iel)*cb26(iel);
    sacasb26  = sa26(iel)*ca26(iel)*sb26(iel);
    
    KK26      = k*[-ca2cb226  -ca2sbcb26 -sacacb26  ca2cb226   ca2sbcb26  sacacb26 ;...
                   -ca2sbcb26 -ca2sb226  -sacasb26  ca2sbcb26  ca2sb226   sacasb26 ;...
                   -sacacb26  -sacasb26  -sa226     sacacb26   sacasb26   sa226    ;...
                    ca2cb226   ca2sbcb26  sacacb26 -ca2cb226  -ca2sbcb26 -sacacb26 ;...
                    ca2sbcb26  ca2sb226   sacasb26 -ca2sbcb26 -ca2sb226  -sacasb26 ;...
                    sacacb26   sacasb26   sa226    -sacacb26  -sacasb26  -sa226   ];  % Stiffness matrix for bar vector 2-6
    
    ca2cb237  = ca37(iel)^2*cb37(iel)^2;
    ca2sb237  = ca37(iel)^2*sb37(iel)^2;
    sa237     = sa37(iel)^2;
    ca2sbcb37 = ca37(iel)^2*sb37(iel)*cb37(iel);
    sacacb37  = sa37(iel)*ca37(iel)*cb37(iel);
    sacasb37  = sa37(iel)*ca37(iel)*sb37(iel);
    
    KK37      = k*[-ca2cb237  -ca2sbcb37 -sacacb37  ca2cb237   ca2sbcb37  sacacb37 ;...
                   -ca2sbcb37 -ca2sb237  -sacasb37  ca2sbcb37  ca2sb237   sacasb37 ;...
                   -sacacb37  -sacasb37  -sa237     sacacb37   sacasb37   sa237    ;...
                    ca2cb237   ca2sbcb37  sacacb37 -ca2cb237  -ca2sbcb37 -sacacb37 ;...
                    ca2sbcb37  ca2sb237   sacasb37 -ca2sbcb37 -ca2sb237  -sacasb37 ;...
                    sacacb37   sacasb37   sa237    -sacacb37  -sacasb37  -sa237   ];  % Stiffness matrix for bar vector 3-7
    
    ca2cb248  = ca48(iel)^2*cb48(iel)^2;
    ca2sb248  = ca48(iel)^2*sb48(iel)^2;
    sa248     = sa48(iel)^2;
    ca2sbcb48 = ca48(iel)^2*sb48(iel)*cb48(iel);
    sacacb48  = sa48(iel)*ca48(iel)*cb48(iel);
    sacasb48  = sa48(iel)*ca48(iel)*sb48(iel);
    
    KK48      = k*[-ca2cb248  -ca2sbcb48 -sacacb48  ca2cb248   ca2sbcb48  sacacb48 ;...
                   -ca2sbcb48 -ca2sb248  -sacasb48  ca2sbcb48  ca2sb248   sacasb48 ;...
                   -sacacb48  -sacasb48  -sa248     sacacb48   sacasb48   sa248    ;...
                    ca2cb248   ca2sbcb48  sacacb48 -ca2cb248  -ca2sbcb48 -sacacb48 ;...
                    ca2sbcb48  ca2sb248   sacasb48 -ca2sbcb48 -ca2sb248  -sacasb48 ;...
                    sacacb48   sacasb48   sa248    -sacacb48  -sacasb48  -sa248   ];  % Stiffness matrix for bar vector 4-8
    
    L_0_CrossBars              = L0_CrossBars(iel);
    dofs_CrossBars([1 4 7 10]) = 3*EL2NOD(iel,:)-2;
    dofs_CrossBars([2 5 8 11]) = 3*EL2NOD(iel,:)-1;
    dofs_CrossBars([3 6 9 12]) = 3*EL2NOD(iel,:);
    rows                       = dofs_CrossBars' * ones(1,12);
    cols                       = ones(12,1) * dofs_CrossBars;
        
    % net force for each element (tetrahedron) due to crossbars is given by:
    %                      _                                                                                    _
    %                     |                                                                                      |
    %                     |       [  cacb15  ]          [  cacb26  ]          [  cacb37  ]          [  cacb48  ] |
    %                     |       |          |          |          |          |          |          |          | |
    %                     |       |  casb15  |          |  casb26  |          |  casb37  |          |  casb48  | |
    %                     |       |          |          |          |          |          |          |          | |
    %                     |       |   sa15   |          |   sa26   |          |   sa37   |          |   sa48   | |
    % {fel_Crossbars} = k*|[TT15]'|          | + [TT26]'|          | + [TT37]'|          | + [TT48]'|          | | * L_0
    %                     |       | -cacb15  |          | -cacb26  |          | -cacb37  |          | -cacb48  | |
    %                     |       |          |          |          |          |          |          |          | |
    %                     |       | -casb15  |          | -casb26  |          | -casb37  |          | -casb48  | |
    %                     |       |          |          |          |          |          |          |          | |
    %                     |       [   -sa15  ]          [   -sa26  ]          [   -sa37  ]          [   -sa48  ] |
    %                     |_                                                                                    _|
    
    fel_CrossBars = k*( TT15'*[ ca15(iel) .* cb15(iel);  ...
                                ca15(iel) .* sb15(iel);  ...
                                             sa15(iel);  ...
                               -ca15(iel) .* cb15(iel);  ...
                               -ca15(iel) .* sb15(iel);  ...
                                            -sa15(iel) ] ...
                       +TT26'*[ ca26(iel) .* cb26(iel);  ...
                                ca26(iel) .* sb26(iel);  ...
                                             sa26(iel);  ...
                               -ca26(iel) .* cb26(iel);  ...
                               -ca26(iel) .* sb26(iel);  ...
                                            -sa26(iel) ] ...
                       +TT37'*[ ca37(iel) .* cb37(iel);  ...
                                ca37(iel) .* sb37(iel);  ...
                                             sa37(iel);  ...
                               -ca37(iel) .* cb37(iel);  ...
                               -ca37(iel) .* sb37(iel);  ...
                                            -sa37(iel) ] ...
                       +TT48'*[ ca48(iel) .* cb48(iel);  ...
                                ca48(iel) .* sb48(iel);  ...
                                             sa48(iel);  ...
                               -ca48(iel) .* cb48(iel);  ...
                               -ca48(iel) .* sb48(iel);  ...
                                            -sa48(iel) ] )*L_0_CrossBars;
    
    % stiffness matrix for each element (tetrahedron) due to crossbars is given by:
    %                     _                                                                                   _
    %                    |                                                                                     |
    % [KKel_Crossbars] = |[TT15]'[KK15][TT15] + [TT26]'[KK26][TT26] + [TT37]'[KK37][TT37] + [TT48]'[KK48][TT48]|
    %                    |_                                                                                   _|
    
    KKel_CrossBars = TT15'*KK15*TT15 + TT26'*KK26*TT26 + TT37'*KK37*TT37 + TT48'*KK48*TT48;
    
    KK1_CrossBars(:,iel) = KKel_CrossBars(:);
    KKi_CrossBars(:,iel) = rows(:);
    KKj_CrossBars(:,iel) = cols(:);
    F_le(dofs_CrossBars) = F_le(dofs_CrossBars) + fel_CrossBars;
end

end % END OF SUBFUNCTION crossbars