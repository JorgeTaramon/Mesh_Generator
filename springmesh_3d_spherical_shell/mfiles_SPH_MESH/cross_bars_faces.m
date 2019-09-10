function [KK1_CrossBars_faces,KKi_CrossBars_faces,KKj_CrossBars_faces,F_le] = ...
    cross_bars_faces(MESH,F_le,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: [KK1_CrossBars_faces,KKi_CrossBars_faces,KKj_CrossBars_faces,F_le] = ...
%   cross_bars_faces(MESH,F_le,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%              3                We add 3 new bars (springs) in every face 
%              .                (i.e., springs 1-5,2-6 and 3-4) of each
%             /.\               tetrahedron, asking them to be the same 
%            / . \              length, therefore these 12 springs will  
%           /  .  \             help to generate a better tetrahedron 
%          /   .   \            (close to a regular one).
%         /    .    \           We are going to generate the 3 new nodes in
%       6.     .     .5         each face: 4,5,6, but we are not going to 
%       /   .  .  .   \         add new dofs. The (x,y,z) for nodes 4,5,6 
%      /       .       \        are only used to calculate the direction 
%     /     .  .  .     \       cosines for the rotation matrix. After 
%    /   .     .     .   \      that, all the information related to 4,5,6 
%   / .        .        . \     will be divided and then added onto nodes 
%  .-----------.-----------.    1,2,3. For example: for spring 1-5, the 
% 1            4            2   information related to node 1 will be fully
%                               added to node 1, but the info related to 
%  node 5 will be splited into 2 halves, and be added to node 2 and node 3.
%
% Input:
%   MESH                : [structure]     : structure containing the mesh
%   F_le                : [column vector] : force vector term due to the 
%                                           effect of a preferred 
%                                           zero-deformation element length
%   GUIDE_MESH          : [structure]     : structure containing guide mesh
%   INTERFACE           : [structure]     : structure containing interface
%                                           settings 
%   SETTINGS            : [structure]     : structure containing mesh
%                                           settings 
% Output:
%   KK1_CrossBars_faces : [matrix]        : numerical value for the element
%                                           (i,j) in the stiffness matrix
%   KKi_CrossBars_faces : [matrix]        : i index (arrow) in the 
%                                           stiffness matrix 
%   KKj_CrossBars_faces : [matrix]        : j index (column) in the 
%                                           stiffness matrix
%   F_le                : [column vector] : force vector term due to the  
%                                           effect of a preferred 
%                                           zero-deformation element length
%                                           after applying cross bars 
% JMT Jan 2016
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

GCOORD             = MESH.GCOORD;
EL2NOD             = MESH.EL2NOD;
k                  = 1; % define spring constant (by default is 1)
faces              = [EL2NOD(:,[1,2,3]); ...
                      EL2NOD(:,[1,2,4]); ...
                      EL2NOD(:,[1,3,4]); ...
                      EL2NOD(:,[2,3,4])];
faces_sorted       = sort(faces,2); % sort nodes of the faces
[faces_unique,~,~] = unique(faces_sorted,'rows');
n_faces            = size(faces_unique,1);
bars_faces         = [faces_unique(:,[1,2]);faces_unique(:,[3,1]);faces_unique(:,[2,3])];
[~,~,J_faces]      = unique(sort(bars_faces,2),'rows'); % J_faces will be used later

bars               = [EL2NOD(:,[1,2]); ...
                      EL2NOD(:,[1,3]); ...
                      EL2NOD(:,[1,4]); ...
                      EL2NOD(:,[2,3]); ...
                      EL2NOD(:,[2,4]); ...
                      EL2NOD(:,[3,4])]; % bars defining each tetrahedron
barmid             = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoint of each bar
% compute L0 for all bars (even those bars whose ends are fixed nodes)
switch SETTINGS.refinement
    case 'regular'
        L0 = SETTINGS.h0*ones(size(barmid,1),1);
    case 'interface'
        L0 = bar_L0_interface(barmid,SETTINGS.h0,INTERFACE);
    case 'guide_mesh'
        L0 = bar_L0_guide(barmid,GUIDE_MESH);
end
% List of bar vectors for each cross bar of the face
barvec15 = (GCOORD(faces_unique(:,2),:) + GCOORD(faces_unique(:,3),:))/2 - GCOORD(faces_unique(:,1),:); 
barvec26 = (GCOORD(faces_unique(:,1),:) + GCOORD(faces_unique(:,3),:))/2 - GCOORD(faces_unique(:,2),:);
barvec34 = (GCOORD(faces_unique(:,1),:) + GCOORD(faces_unique(:,2),:))/2 - GCOORD(faces_unique(:,3),:);
% Generate all the sin(alpha),cos(alpha),sin(betha),cos(betha) for each bar
alpha15 = atan2d(barvec15(:,3),sqrt(barvec15(:,1).^2+barvec15(:,2).^2));
beta15  = atan2d(barvec15(:,2),barvec15(:,1));
sa15    = sind(alpha15);
ca15    = cosd(alpha15);
sb15    = sind(beta15);
cb15    = cosd(beta15);

alpha26 = atan2d(barvec26(:,3),sqrt(barvec26(:,1).^2+barvec26(:,2).^2));
beta26  = atan2d(barvec26(:,2),barvec26(:,1));
sa26    = sind(alpha26);
ca26    = cosd(alpha26);
sb26    = sind(beta26);
cb26    = cosd(beta26);

alpha34 = atan2d(barvec34(:,3),sqrt(barvec34(:,1).^2+barvec34(:,2).^2));
beta34  = atan2d(barvec34(:,2),barvec34(:,1));
sa34    = sind(alpha34);
ca34    = cosd(alpha34);
sb34    = sind(beta34);
cb34    = cosd(beta34);

% Generate all the L0 info: the 3 cross bars in each face share a same L0, and this L0 = sqrt(3)/2*average(L0 value for the edges)
L0_CrossBars_faces = sqrt(3)/2*mean(reshape(L0(J_faces),n_faces,3),2); % n_faces-by-1
%---------------------------- what the above line does -------------------------------
% L0_with_duplication   = MESH.L0(J);
% L0_reshape            = reshape(L0_with_duplication,n_faces,3); % n_faces-by-3 matrix
% L0_mean_for_each_face = mean(L0_reshape,2);                     % n_faces-by-1 vector
% L0_CrossBars_faces    = sqrt(3)/2*L0_mean_for_each_face;
%-------------------------------------------------------------------------------------

KK1_CrossBars_faces  = zeros(81,n_faces);
KKi_CrossBars_faces  = zeros(81,n_faces);
KKj_CrossBars_faces  = zeros(81,n_faces);
dofs_CrossBars_faces = zeros(1,9);

TT15 =  [1  0  0  0  0  0  0  0  0  ;...
         0  1  0  0  0  0  0  0  0  ;...
         0  0  1  0  0  0  0  0  0  ;...
         0  0  0 1/2 0  0 1/2 0  0  ;...
         0  0  0  0 1/2 0  0 1/2 0  ;...
         0  0  0  0  0 1/2 0  0 1/2 ];

TT26 =  [0  0  0  1  0  0  0  0  0  ;...
         0  0  0  0  1  0  0  0  0  ;...
         0  0  0  0  0  1  0  0  0  ;...
        1/2 0  0  0  0  0 1/2 0  0  ;...
         0 1/2 0  0  0  0  0 1/2 0  ;...
         0  0 1/2 0  0  0  0  0 1/2 ];

TT34 =  [0  0  0  0  0  0  1  0  0  ;...
         0  0  0  0  0  0  0  1  0  ;...
         0  0  0  0  0  0  0  0  1  ;...
        1/2 0  0 1/2 0  0  0  0  0  ;...
         0 1/2 0  0 1/2 0  0  0  0  ;...
         0  0 1/2 0  0 1/2 0  0  0  ];

% loop over all the faces, calculate the fel and KK info, also generate the
% indices for all the fel, KK values (for sparse matrix)
for iface = 1:n_faces
    
    ca2cb2_15  = ca15(iface)^2*cb15(iface)^2;
    ca2sb2_15  = ca15(iface)^2*sb15(iface)^2;
    sa2_15     = sa15(iface)^2;
    ca2sbcb_15 = ca15(iface)^2*sb15(iface)*cb15(iface);
    sacacb_15  = sa15(iface)*ca15(iface)*cb15(iface);
    sacasb_15  = sa15(iface)*ca15(iface)*sb15(iface);
    
    KK15       = k*[-ca2cb2_15  -ca2sbcb_15 -sacacb_15  ca2cb2_15   ca2sbcb_15  sacacb_15 ;...
                    -ca2sbcb_15 -ca2sb2_15  -sacasb_15  ca2sbcb_15  ca2sb2_15   sacasb_15 ;...
                    -sacacb_15  -sacasb_15  -sa2_15     sacacb_15   sacasb_15   sa2_15    ;...
                     ca2cb2_15   ca2sbcb_15  sacacb_15 -ca2cb2_15  -ca2sbcb_15 -sacacb_15 ;...
                     ca2sbcb_15  ca2sb2_15   sacasb_15 -ca2sbcb_15 -ca2sb2_15  -sacasb_15 ;...
                     sacacb_15   sacasb_15   sa2_15    -sacacb_15  -sacasb_15  -sa2_15   ];   % Stiffness matrix for bar vector 1-5
    
    ca2cb2_26  = ca26(iface)^2*cb26(iface)^2;
    ca2sb2_26  = ca26(iface)^2*sb26(iface)^2;
    sa2_26     = sa26(iface)^2;
    ca2sbcb_26 = ca26(iface)^2*sb26(iface)*cb26(iface);
    sacacb_26  = sa26(iface)*ca26(iface)*cb26(iface);
    sacasb_26  = sa26(iface)*ca26(iface)*sb26(iface);
    
    KK26       = k*[-ca2cb2_26  -ca2sbcb_26 -sacacb_26  ca2cb2_26   ca2sbcb_26  sacacb_26 ;...
                    -ca2sbcb_26 -ca2sb2_26  -sacasb_26  ca2sbcb_26  ca2sb2_26   sacasb_26 ;...
                    -sacacb_26  -sacasb_26  -sa2_26     sacacb_26   sacasb_26   sa2_26    ;...
                     ca2cb2_26   ca2sbcb_26  sacacb_26 -ca2cb2_26  -ca2sbcb_26 -sacacb_26 ;...
                     ca2sbcb_26  ca2sb2_26   sacasb_26 -ca2sbcb_26 -ca2sb2_26  -sacasb_26 ;...
                     sacacb_26   sacasb_26   sa2_26    -sacacb_26  -sacasb_26  -sa2_26   ];   % Stiffness matrix for bar vector 2-6
    
    ca2cb2_34  = ca34(iface)^2*cb34(iface)^2;
    ca2sb2_34  = ca34(iface)^2*sb34(iface)^2;
    sa2_34     = sa34(iface)^2;
    ca2sbcb_34 = ca34(iface)^2*sb34(iface)*cb34(iface);
    sacacb_34  = sa34(iface)*ca34(iface)*cb34(iface);
    sacasb_34  = sa34(iface)*ca34(iface)*sb34(iface);
    
    KK34       = k*[-ca2cb2_34  -ca2sbcb_34 -sacacb_34  ca2cb2_34   ca2sbcb_34  sacacb_34 ;...
                    -ca2sbcb_34 -ca2sb2_34  -sacasb_34  ca2sbcb_34  ca2sb2_34   sacasb_34 ;...
                    -sacacb_34  -sacasb_34  -sa2_34     sacacb_34   sacasb_34   sa2_34    ;...
                     ca2cb2_34   ca2sbcb_34  sacacb_34 -ca2cb2_34  -ca2sbcb_34 -sacacb_34 ;...
                     ca2sbcb_34  ca2sb2_34   sacasb_34 -ca2sbcb_34 -ca2sb2_34  -sacasb_34 ;...
                     sacacb_34   sacasb_34   sa2_34    -sacacb_34  -sacasb_34  -sa2_34   ];   % Stiffness matrix for bar vector 3-4
    
    L_0_CrossBar_faces            = L0_CrossBars_faces(iface);
    dofs_CrossBars_faces([1 4 7]) = 3*faces_unique(iface,:)-2; % x-dofs
    dofs_CrossBars_faces([2 5 8]) = 3*faces_unique(iface,:)-1; % y-dofs
    dofs_CrossBars_faces([3 6 9]) = 3*faces_unique(iface,:);   % z-dofs
    rows                          = dofs_CrossBars_faces' * ones(1,9);
    cols                          = ones(9,1) * dofs_CrossBars_faces;
    
    % net force for each face due to crossbars is given by:
    %                            _                                                              _
    %                           |                                                                |
    %                           |       [  cacb15  ]          [  cacb26  ]          [  cacb34  ] |
    %                           |       |          |          |          |          |          | |
    %                           |       |  casb15  |          |  casb26  |          |  casb34  | |
    %                           |       |          |          |          |          |          | |
    %                           |       |    sa15  |          |    sa26  |          |    sa34  | |
    % {fel_CrossBars_faces} = k*|[TT15]'|          | + [TT26]'|          | + [TT34]'|          | | * L_0
    %                           |       | -cacb15  |          | -cacb26  |          | -cacb34  | |
    %                           |       |          |          |          |          |          | |
    %                           |       | -casb15  |          | -casb26  |          | -casb34  | |
    %                           |       |          |          |          |          |          | |
    %                           |       [   -sa15  ]          [   -sa26  ]          [   -sa34  ] |
    %                           |_                                                              _|
    
    fel_CrossBars_faces = k*( TT15'*[ ca15(iface) .* cb15(iface);  ...
                                      ca15(iface) .* sb15(iface);  ...
                                                     sa15(iface);  ...
                                     -ca15(iface) .* cb15(iface);  ...
                                     -ca15(iface) .* sb15(iface);  ...
                                                    -sa15(iface) ] ...
                             +TT26'*[ ca26(iface) .* cb26(iface);  ...
                                      ca26(iface) .* sb26(iface);  ...
                                                     sa26(iface);  ...
                                     -ca26(iface) .* cb26(iface);  ...
                                     -ca26(iface) .* sb26(iface);  ...
                                                    -sa26(iface) ] ...
                             +TT34'*[ ca34(iface) .* cb34(iface);  ...
                                      ca34(iface) .* sb34(iface);  ...
                                                     sa34(iface);  ...
                                     -ca34(iface) .* cb34(iface);  ...
                                     -ca34(iface) .* sb34(iface);  ...
                                                    -sa34(iface) ] )*L_0_CrossBar_faces;
    
    % stiffness matrix for each face due to crossbars is given by:
    %
    % [KKel_CrossBars_faces] = [TT15]'[KK15][TT15] + [TT26]'[KK26][TT26] + [TT34]'[KK34][TT34]
    
    KKel_CrossBars_faces = TT15'*KK15*TT15 + TT26'*KK26*TT26 + TT34'*KK34*TT34;
    
    KK1_CrossBars_faces(:,iface) = KKel_CrossBars_faces(:);
    KKi_CrossBars_faces(:,iface) = rows(:);
    KKj_CrossBars_faces(:,iface) = cols(:);
    F_le(dofs_CrossBars_faces)   = F_le(dofs_CrossBars_faces) + fel_CrossBars_faces;
    
end
end % END OF SUBFUNCTION cross_bars_faces