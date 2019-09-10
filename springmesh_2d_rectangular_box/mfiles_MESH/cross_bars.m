function [KK1_CrossBars,KKi_CrossBars,KKj_CrossBars,F_le] = cross_bars(MESH,F_le)
% Usage: [KK1_CrossBars,KKi_CrossBars,KKj_CrossBars,F_le] = cross_bars(MESH,F_le)
%
% Purpose:
%
%              3                We add 3 new bars (springs) in every  
%              .                triangle (i.e., springs 1-5,2-6 and 3-4) 
%             /.\               asking them to be the same length, 
%            / . \              therefore these 3 springs will help to 
%           /  .  \             generate better triangles (close to 
%          /   .   \            equilateral ones).
%         /    .    \           We are going to generate the 3 new nodes: 
%       6.     .     .5         4,5,6, but we are not going to add new 
%       /   .  .  .   \         dofs. The (x,y) for nodes 4,5,6 are only 
%      /       .       \        used to calculate the direction cosines for
%     /     .  .  .     \       the rotation matrix. After that, all the 
%    /   .     .     .   \      information related to 4,5,6 will be 
%   / .        .        . \     divided and then added onto nodes 1,2,3.
%  .-----------.-----------.    For example: for spring 1-5, the 
% 1            4            2   information related to node 1 will be fully
%                               added to node 1, but the info related to 
%  node 5 will be splited into 2 halves, and be added to node 2 and node 3.
%
% Input:
%   MESH          : [structure]     : structure containing the mesh
%   F_le          : [column vector] : force vector term due to the effect
%                                     of a preferred zero-deformation
%                                     element length
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
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL7
%--------------------------------------------------------------------------

GCOORD   = MESH.GCOORD;
EL2NOD   = MESH.EL2NOD;
k        = 1; % define spring constant (default is 1)
nel_tri  = length(EL2NOD);
bars     = [EL2NOD(:,[1,2]);EL2NOD(:,[2,3]);EL2NOD(:,[3,1])]; % bars defining each triangle
[~,~,J]  = unique(sort(bars,2),'rows');
% List of bar vector for each cross bar of the triangle
barvec15 = (GCOORD(EL2NOD(:,2),:) + GCOORD(EL2NOD(:,3),:))/2 - GCOORD(EL2NOD(:,1),:); % vector from node 1 to node 5
barvec26 = (GCOORD(EL2NOD(:,1),:) + GCOORD(EL2NOD(:,3),:))/2 - GCOORD(EL2NOD(:,2),:); % vector from node 2 to node 6
barvec34 = (GCOORD(EL2NOD(:,1),:) + GCOORD(EL2NOD(:,2),:))/2 - GCOORD(EL2NOD(:,3),:); % vector from node 3 to node 4
% Generate all the sin(alpha),cos(alpha),sin(betha),cos(betha) for each bar
L15      = sqrt(sum(barvec15.^2,2)); % L = bar lengths
L26      = sqrt(sum(barvec26.^2,2));
L34      = sqrt(sum(barvec34.^2,2));
s15      = barvec15(:,2) ./ L15(:); % sin(alpha(bar)) for RR matrix
c15      = barvec15(:,1) ./ L15(:); % cos(alpha(bar)) for RR matrix
s26      = barvec26(:,2) ./ L26(:); % sin(alpha(bar)) for RR matrix
c26      = barvec26(:,1) ./ L26(:); % cos(alpha(bar)) for RR matrix
s34      = barvec34(:,2) ./ L34(:); % sin(alpha(bar)) for RR matrix
c34      = barvec34(:,1) ./ L34(:); % cos(alpha(bar)) for RR matrix

% Generate all the L0 info: the 3 cross bars in each triangle share a same 
% L0, and this L0 = sqrt(3)/2*average(3 L0 value for the edges)
L0_CrossBars = sqrt(3)/2*mean(reshape(MESH.L0(J),nel_tri,3),2); % nel_tri-by-1
%----------------------------- what the above line does ------------------------------
% L0_with_duplication  = MESH.L0(J);
% L0_reshape           = reshape(L0_with_duplication,nel_tri,3); % nel_tri-by-3 matrix
% L0_mean_for_each_tri = mean(L0_reshape,2);                     % nel_tri-by-1 vector
% L0_CrossBars         = sqrt(3)/2*L0_mean_for_each_tri;
%-------------------------------------------------------------------------------------

KK1_CrossBars  = zeros(36,nel_tri);
KKi_CrossBars  = zeros(36,nel_tri);
KKj_CrossBars  = zeros(36,nel_tri);
dofs_CrossBars = zeros(1,6);

TT15 = [1    0    0    0    0    0   ; ...
        0    1    0    0    0    0   ; ...
        0    0    0.5  0    0.5  0   ; ...
        0    0    0    0.5  0    0.5];

TT26 = [0    0    1    0    0    0   ; ...
        0    0    0    1    0    0   ; ...
        0.5  0    0    0    0.5  0   ; ...
        0    0.5  0    0    0    0.5];

TT34 = [0    0    0    0    1    0   ; ...
        0    0    0    0    0    1   ; ...
        0.5  0    0.5  0    0    0   ; ...
        0    0.5  0    0.5  0    0  ];

% loop over all the triangles, calculate the fel and KK info, also 
% generates the indices for all the fel, KK values (for sparse matrix)
for iel = 1:nel_tri
    
    s215 = s15(iel)^2;
    c215 = c15(iel)^2;
    s226 = s26(iel)^2;
    c226 = c26(iel)^2;
    s234 = s34(iel)^2;
    c234 = c34(iel)^2;
    sc15 = s15(iel)*c15(iel);
    sc26 = s26(iel)*c26(iel);
    sc34 = s34(iel)*c34(iel);
    
    KK15 = k*[-c215  -sc15    c215    sc15; ...
              -sc15  -s215    sc15    s215; ...
               c215   sc15   -c215   -sc15; ...
               sc15   s215   -sc15   -s215]; % Stiffness matrix for bar vector 1-5
    
    KK26 = k*[-c226  -sc26    c226    sc26; ...
              -sc26  -s226    sc26    s226; ...
               c226   sc26   -c226   -sc26; ...
               sc26   s226   -sc26   -s226]; % Stiffness matrix for bar vector 2-6
    
    KK34 = k*[-c234  -sc34    c234    sc34; ...
              -sc34  -s234    sc34    s234; ...
               c234   sc34   -c234   -sc34; ...
               sc34   s234   -sc34   -s234]; % Stiffness matrix for bar vector 3-4
    
    L_0_CrossBars           = L0_CrossBars(iel);
    dofs_CrossBars([1 3 5]) = 2*EL2NOD(iel,:)-1; % x-dofs
    dofs_CrossBars([2 4 6]) = 2*EL2NOD(iel,:);   % y-dofs
    rows                    = dofs_CrossBars' * ones(1,6);
    cols                    = ones(6,1) * dofs_CrossBars;
        
    % net force for each element (triangle) due to crossbars is given by:
    %                      _                                                  _
    %                     |                                                    |
    %                     |       [  c15 ]          [  c26 ]          [  c34 ] |
    %                     |       |      |          |      |          |      | |
    %                     |       |  s15 |          |  s26 |          |  s34 | |
    % {fel_Crossbars} = k*|[TT15]'|      | + [TT26]'|      | + [TT34]'|      | | * L_0
    %                     |       | -c15 |          | -c26 |          | -c34 | |
    %                     |       |      |          |      |          |      | |
    %                     |       [ -s15 ]          [ -s26 ]          [ -s34 ] |
    %                     |_                                                  _|
    
    fel_CrossBars = k*(TT15'*[c15(iel);s15(iel);-c15(iel);-s15(iel)]...
                      +TT26'*[c26(iel);s26(iel);-c26(iel);-s26(iel)]...
                      +TT34'*[c34(iel);s34(iel);-c34(iel);-s34(iel)])*L_0_CrossBars;
    
    % stiffness matrix for each element (triangle) due to crossbars is given by:
    %                     _                                                             _
    %                    |                                                               |
    % [KKel_Crossbars] = |[TT15]'[KK15][TT15] + [TT26]'[KK26][TT26] + [TT34]'[KK34][TT34]|
    %                    |_                                                             _|
    
    KKel_CrossBars = TT15'*KK15*TT15+TT26'*KK26*TT26+TT34'*KK34*TT34;
    
    KK1_CrossBars(:,iel) = KKel_CrossBars(:);
    KKi_CrossBars(:,iel) = rows(:);
    KKj_CrossBars(:,iel) = cols(:);
    F_le(dofs_CrossBars) = F_le(dofs_CrossBars) + fel_CrossBars;
end

end % END OF SUBFUNCTION crossbars