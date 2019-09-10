function [x_I,y_I,z_I,r_I,n1,n2,n3,n4] = interior_sphere(GCOORD,EL2NOD,SETTINGS)
% Usage: [x_I,y_I,z_I,r_I,n1,n2,n3,n4] = interior_sphere(GCOORD,EL2NOD,SETTINGS)
%
% Purpose:
%   For any tetrahedron compute the coordinates of the centre of the
%   insphere (x_I,y_I,z_I), the radius (r_I) and the unit normal vectors to
%   each face (n1,n2,n3,n4).
%   Condition that defines the insphere: The perpendicular distance from
%   the incenter to the faces of the tetrahedron are all equal.
%   For a given face, f, let n its unit normal vector pointing out of the
%   tetrahedron. Let (x_I,y_y,z_I) be the unknown center of the insphere 
%   and r_I its radius. We look for a condition that (x_I,y_y,z_I) + r_I*n 
%   lies on the face f. Let v12 and v13 be the displacement vectors of two 
%   edges of f, and let P0 be one of the vertices of f. Then the unit  
%   outward normal n will be the normalization of ±(v12 x v13). The standar 
%   point-normal form for the equation of a plane containing a point P0 and
%   perpendicular to a vector n is, for an arbitrary point P on the plane:
%               
%                              (P-P0)*n = 0
%
%   which is equivalent to
%
%               [(x_I,y_y,z_I) + r_I*n_i - (x_i,y_i,z_i)]*n_i
% 
%   where i = 1 to 4 (tetrahedron faces).
%   To get the orientation of the normal n right, we need to compute for
%   each face the sign of the minus dot product between the cross product
%   of the vectors that define the face and the third vector which is going
%   towards the interior of the face, e,g for the face 1:
% 
%                          s1 = -(v12 x v13)*v14
%
%   So the unit normal vetor pointing outward the tetrahedron is given by:
%
%                 n1 = sign(s1)*(v12 x v13)/||(v12 x v13)||
%
%   In general for each tetrahedron:
%
%                   4                           
%                   .                          
%                  / \\                        
%                 /   \  \                     
%                /     \    \                  
%               /       \      \               
%              /         \        \            
%             /           \          \         
%            /             \            \      
%           /               \       _ _ _. 3 
%          /                _\_ _---   /       
%         /        _ _ _----  \      /         
%        / _ _ _---            \   /          
%       .-----------------------./           
%      1                         2           
%
%
%          |            |  vectors |      scale factor     |                             | 
%    Face  |     P0     | defining |     for orientation   |      unit normal vector     | 
%          |            | the face |    of normal vector   |                             | 
%   -------|------------|----------|-----------------------|-----------------------------| 
%          |            |          |                       |       sign(s1)*(v12 x v13)  | 
%      1   | (x1,y1,z1) | v12, v13 | s1 = -(v12 x v13)*v14 | n1 = ---------------------- | 
%          |            |          |                       |         ||(v12 x v13)||     | 
%   -------|------------|----------|-----------------------|-----------------------------| 
%          |            |          |                       |       sign(s2)*(v23 x v24)  | 
%      2   | (x2,y2,z2) | v23, v24 | s2 = -(v23 x v24)*v21 | n2 = ---------------------- | 
%          |            |          |                       |         ||(v23 x v24)||     | 
%   -------|------------|----------|-----------------------|-----------------------------| 
%          |            |          |                       |       sign(s3)*(v31 x v34)  | 
%      3   | (x3,y3,z3) | v31, v34 | s3 = -(v31 x v34)*v32 | n3 = ---------------------- | 
%          |            |          |                       |         ||(v31 x v34)||     | 
%   -------|------------|----------|-----------------------|-----------------------------| 
%          |            |          |                       |       sign(s4)*(v41 x v42)  | 
%      4   | (x4,y4,z4) | v41, v42 | s4 = -(v41 x v42)*v43 | n4 = ---------------------- | 
%          |            |          |                       |         ||(v41 x v42)||     | 
%   -------|------------|----------|-----------------------|-----------------------------| 
%
% Input:
%   GCOORD : [matrix] : coorfinates of the nodes of the mesh
%   EL2NOD : [matrix] : conectivity matrix of the mesh
%
% Output:
%   x_I    : [vector] : x-coordinate of incentre
%   y_I    : [vector] : y-coordinate of incentre
%   z_I    : [vector] : z-coordinate of incentre
%   r_I    : [vector] : radius of the insphere
%   n1     : [vector] : unit normal vector to the face 1 pointing outward
%                       the tetrahedron 
%   n2     : [vector] : unit normal vector to the face 2 pointing outward
%                       the tetrahedron 
%   n3     : [vector] : unit normal vector to the face 3 pointing outward
%                       the tetrahedron 
%   n4     : [vector] : unit normal vector to the face 4 pointing outward
%                       the tetrahedron 
%                 
% JMT Oct 2015
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% COMPUTE VECTORS THAT DEFINE THE FACES AND ITS UNIT NORMAL VECTORS 
%==========================================================================
nel        = size(EL2NOD,1); % number of tetrahedrons
bars       = [EL2NOD(:,[1,2]); ...
              EL2NOD(:,[1,3]); ...
              EL2NOD(:,[1,4]); ...
              EL2NOD(:,[2,3]); ...
              EL2NOD(:,[2,4]); ...
              EL2NOD(:,[3,4])]; % bars defining each tetrahedron
barvec     = GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:); % vectors v12,v13,v14,v23,v24,v34
barvec_inv = GCOORD(bars(:,1),:) - GCOORD(bars(:,2),:); % vectors v21,v31,v41,v32,v42,v43                                  

v12 = barvec(1:nel,:);
v13 = barvec(nel+1:2*nel,:);
v14 = barvec(2*nel+1:3*nel,:);
v23 = barvec(3*nel+1:4*nel,:);
v24 = barvec(4*nel+1:5*nel,:);
v34 = barvec(5*nel+1:6*nel,:);

v21 = barvec_inv(1:nel,:);
v31 = barvec_inv(nel+1:2*nel,:);
v41 = barvec_inv(2*nel+1:3*nel,:);
v32 = barvec_inv(3*nel+1:4*nel,:);
v42 = barvec_inv(4*nel+1:5*nel,:);
v43 = barvec_inv(5*nel+1:6*nel,:);

s1  = sign(-dot(cross(v12,v13,2),v14,2)); % factor for orientation of unit normal vector in face 1 
s2  = sign(-dot(cross(v23,v24,2),v21,2)); % factor for orientation of unit normal vector in face 2 
s3  = sign(-dot(cross(v31,v34,2),v32,2)); % factor for orientation of unit normal vector in face 3 
s4  = sign(-dot(cross(v41,v42,2),v43,2)); % factor for orientation of unit normal vector in face 4 

n1  = repmat(s1,1,3).*cross(v12,v13,2)./repmat(sqrt(sum(cross(v12,v13,2).^2,2)),1,3); % unit normal vector in face 1 pointing outward the tetrahedron
n2  = repmat(s2,1,3).*cross(v23,v24,2)./repmat(sqrt(sum(cross(v23,v24,2).^2,2)),1,3); % unit normal vector in face 2 pointing outward the tetrahedron
n3  = repmat(s3,1,3).*cross(v14,v13,2)./repmat(sqrt(sum(cross(v14,v13,2).^2,2)),1,3); % unit normal vector in face 3 pointing outward the tetrahedron
n4  = repmat(s4,1,3).*cross(v12,v14,2)./repmat(sqrt(sum(cross(v12,v14,2).^2,2)),1,3); % unit normal vector in face 4 pointing outward the tetrahedron

a1  = n1(:,1);
a2  = n2(:,1);
a3  = n3(:,1);
a4  = n4(:,1);

b1  = n1(:,2);
b2  = n2(:,2);
b3  = n3(:,2);
b4  = n4(:,2);

c1  = n1(:,3);
c2  = n2(:,3);
c3  = n3(:,3);
c4  = n4(:,3);

d1  = dot(GCOORD(EL2NOD(:,1),:),n1,2);
d2  = dot(GCOORD(EL2NOD(:,2),:),n2,2);
d3  = dot(GCOORD(EL2NOD(:,3),:),n3,2);
d4  = dot(GCOORD(EL2NOD(:,4),:),n4,2);

%==========================================================================
% SOLVE FOR THE INCENTRE AND INRADIUS
%==========================================================================
switch SETTINGS.solver_balloon_forces
    case 'gauss'
        % Solving by analytical solution using Gauss Method
        z_I = ((((a2-a1).*(b3-b1)-(a3-a1).*(b2-b1)).*((a2-a1).*(d4-d1)-(a4-a1).*(d2-d1)))  -  ...
               (((a2-a1).*(b4-b1)-(a4-a1).*(b2-b1)).*((a2-a1).*(d3-d1)-(a3-a1).*(d2-d1)))) ./ ...
              ((((a2-a1).*(b3-b1)-(a3-a1).*(b2-b1)).*((a2-a1).*(c4-c1)-(a4-a1).*(c2-c1)))  -  ...
               (((a2-a1).*(b4-b1)-(a4-a1).*(b2-b1)).*((a2-a1).*(c3-c1)-(a3-a1).*(c2-c1))));
        
        y_I =       ((a2-a1).*(d3-d1)-(a3-a1).*(d2-d1)) ./ ...
                    ((a2-a1).*(b3-b1)-(a3-a1).*(b2-b1)) -  ...
               z_I.*((a2-a1).*(c3-c1)-(a3-a1).*(c2-c1)) ./ ...
                    ((a2-a1).*(b3-b1)-(a3-a1).*(b2-b1));
        x_I = (d2-d1)./(a2-a1) - z_I.*(c2-c1)./(a2-a1) - y_I.*(b2-b1)./(a2-a1);
        r_I = d1 - z_I.*c1 - y_I.*b1 - x_I.*a1;
        
    case '3D_matrix'
        % Solving by 3-D matrix with "for" loop
        M        = zeros(4,4,nel);
        M(1,1,:) = 1;
        M(1,2,:) = a1;
        M(1,3,:) = b1;
        M(1,4,:) = c1;
        M(2,1,:) = 1;
        M(2,2,:) = a2;
        M(2,3,:) = b2;
        M(2,4,:) = c2;
        M(3,1,:) = 1;
        M(3,2,:) = a3;
        M(3,3,:) = b3;
        M(3,4,:) = c3;
        M(4,1,:) = 1;
        M(4,2,:) = a4;
        M(4,3,:) = b4;
        M(4,4,:) = c4;
        
        rhs        = zeros(4,1,nel);
        rhs(1,1,:) = d1;
        rhs(2,1,:) = d2;
        rhs(3,1,:) = d3;
        rhs(4,1,:) = d4;
        
        x = zeros(4*nel,1);
        
        for i = 1:nel
            x(4*i-3:4*i) = M(:,:,i)\rhs(:,:,i);  % unknows: coordinates of incentre (x_I,y_I,z_I) and inradius (r_I) for each element
        end
        
        r_I = x(1:4:end);
        x_I = x(2:4:end);
        y_I = x(3:4:end);
        z_I = x(4:4:end);
        
    case 'sparse'
        % Solving by sparse matrix
        % Pre-allocate memory: number of non-zero slots
        nnz_M          = 16*nel; % 4*nel + 12*nel
        KKi_M          = zeros(nnz_M,1);
        KKj_M          = KKi_M;
        KK1_M          = ones(nnz_M,1);
        KKi_M(1:4*nel) = (1:4*nel)';
        KKj_M(1:4*nel) = (1:4*nel)';
        
        % Diagonal slots: a1,b2,c3,d4 for each element
        
        KK1_M(1:4:4*nel-3) = 1;  % elements (1,1),(5,5),(9,9)...
        KK1_M(2:4:4*nel-2) = a2; % elements (2,2),(6,6),(10,10)...
        KK1_M(3:4:4*nel-1) = b3; % elements (3,3),(7,7),(11,11)...
        KK1_M(4:4:4*nel)   = c4; % elements (3,3),(8,8),(12,12)...
        
        % Off-diagonal slots: b1,c1,d1,a2,c2,d2,a3,b3,d3,a4,b4,c4 for each element
        
        KKi_M(4*nel+1:5*nel) = (1:4:4*nel-3);
        KKj_M(4*nel+1:5*nel) = (2:4:4*nel-2);
        KK1_M(4*nel+1:5*nel) = a1; % elements (1,2),(5,6),(9,10)...
        
        KKi_M(5*nel+1:6*nel) = (1:4:4*nel-3);
        KKj_M(5*nel+1:6*nel) = (3:4:4*nel-1);
        KK1_M(5*nel+1:6*nel) = b1; % elements (1,3),(5,7),(9,11)...
        
        KKi_M(6*nel+1:7*nel) = (1:4:4*nel-3);
        KKj_M(6*nel+1:7*nel) = (4:4:4*nel  );
        KK1_M(6*nel+1:7*nel) = c1; % elements (1,4),(5,8),(9,12)...
        
        KKi_M(7*nel+1:8*nel) = (2:4:4*nel-2);
        KKj_M(7*nel+1:8*nel) = (1:4:4*nel-3);
        KK1_M(7*nel+1:8*nel) = 1; % elements (2,1),(6,5),(10,9)...
        
        KKi_M(8*nel+1:9*nel) = (2:4:4*nel-2);
        KKj_M(8*nel+1:9*nel) = (3:4:4*nel-1);
        KK1_M(8*nel+1:9*nel) = b2; % elements (2,3),(6,7),(10,11)...
        
        KKi_M(9*nel+1:10*nel) = (2:4:4*nel-2);
        KKj_M(9*nel+1:10*nel) = (4:4:4*nel  );
        KK1_M(9*nel+1:10*nel) = c2; % elements (2,4),(6,8),(10,12)...
        
        KKi_M(10*nel+1:11*nel) = (3:4:4*nel-1);
        KKj_M(10*nel+1:11*nel) = (1:4:4*nel-3);
        KK1_M(10*nel+1:11*nel) = 1; % elements (3,1),(7,5),(11,9)...
        
        KKi_M(11*nel+1:12*nel) = (3:4:4*nel-1);
        KKj_M(11*nel+1:12*nel) = (2:4:4*nel-2);
        KK1_M(11*nel+1:12*nel) = a3; % elements (3,2),(7,6),(11,10)...
        
        KKi_M(12*nel+1:13*nel) = (3:4:4*nel-1);
        KKj_M(12*nel+1:13*nel) = (4:4:4*nel  );
        KK1_M(12*nel+1:13*nel) = c3; % elements (3,4),(7,8),(11,12)...
        
        KKi_M(13*nel+1:14*nel) = (4:4:4*nel  );
        KKj_M(13*nel+1:14*nel) = (1:4:4*nel-3);
        KK1_M(13*nel+1:14*nel) = 1; % elements (4,1),(8,5),(12,9)...
        
        KKi_M(14*nel+1:15*nel) = (4:4:4*nel  );
        KKj_M(14*nel+1:15*nel) = (2:4:4*nel-2);
        KK1_M(14*nel+1:15*nel) = a4; % elements (4,2),(8,6),(12,10)...
        
        KKi_M(15*nel+1:16*nel) = (4:4:4*nel  );
        KKj_M(15*nel+1:16*nel) = (3:4:4*nel-1);
        KK1_M(15*nel+1:16*nel) = b4; % elements (4,3),(8,7),(12,11)...
        
        M   = sparse2(KKi_M(:),KKj_M(:),KK1_M(:));
        X   = zeros(4*nel,1); % unknows: coordinates of incentre (x_I,y_I,z_I) and inradius (r_I) for each element
        rhs = reshape([d1 d2 d3 d4]',4*nel,1);
        
        X(:) = M\rhs;
        
        r_I = X(1:4:end);
        x_I = X(2:4:end);
        y_I = X(3:4:end);
        z_I = X(4:4:end);
end
end % END OF SUBFUNCTION interior_sphere