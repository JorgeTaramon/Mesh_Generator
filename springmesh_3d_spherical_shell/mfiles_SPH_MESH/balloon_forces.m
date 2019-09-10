function F_balloon = balloon_forces(MESH,EL2NOD_sorted,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: F_balloon = balloon_forces(MESH,EL2NOD_sorted,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   This routine is applied in those elements which quality factor is
%   smaller than the quality tolerance (q < q_balloon).
%   Compute the needed force at each node (for each 'bad' element) taking
%   into account the difference between the insphere for the desired
%   tetrahedron and the insphere for actual tetrahedron. This force will 
%   make better elements.
%
% Input:
%   MESH          : [structure]     : structure containing the mesh
%   EL2NOD_sorted : [matrix]        : conectivity matrix for those elements 
%                                     which q < q_balloon 
%   GUIDE_MESH    : [structure]     : structure containing guide mesh
%   INTERFACE     : [structure]     : structure containing interface
%                                     settings 
%   SETTINGS      : [structure]     : structure containing mesh settings
%
% Output:
%   F_balloon     : [column vector] : forces acting at each node of the
%                                     'bad' elements
% JMT Oct 2015
% JMT May 2016: cleaned up
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if nargin == 0
    MESH.GCOORD        = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1; 0 0 0.2];
    MESH.EL2NOD        = delaunay(MESH.GCOORD);
    q_balloon          = 0.6;
    q                  = tetra_mesh_quality(MESH.GCOORD,MESH.EL2NOD);
    [q_sorted,I]       = sort(q);                                % sort q from worst q to best q
    EL2NOD_sorted      = MESH.EL2NOD(I(q_sorted < q_balloon),:); % select those elements which q < q_balloon
    SETTINGS.show_figs = 1;
    SETTINGS.solver_balloon_forces = 'sparse';
end

%==========================================================================
% COMPUTE THE INCENTER I(x_I,y_I,z_I), THE INRADIUS (ACTUAL,r_I, AND
% DESIRED,r0_mean) AND THE FORCES ACTING AT EACH TANGENT POINT OF EACH
% ELEMENT
%==========================================================================
GCOORD      = MESH.GCOORD;
nel         = size(EL2NOD_sorted,1); % number of tetrahedrons
bars        = [EL2NOD_sorted(:,[1,2]); ...
               EL2NOD_sorted(:,[1,3]); ...
               EL2NOD_sorted(:,[1,4]); ...
               EL2NOD_sorted(:,[2,3]); ...
               EL2NOD_sorted(:,[2,4]); ...
               EL2NOD_sorted(:,[3,4])]; % bars defining each tetrahedron
barmid      = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoint of each bar element
[x_I,y_I,z_I,r_I,n1,n2,n3,n4] = interior_sphere(GCOORD,EL2NOD_sorted,SETTINGS);
if nargin == 0
    L0      = [3 3 3 3 3 3; ...
               3 3 3 3 3 3; ...
               3 3 3 3 3 3; ...
               3 3 3 3 3 3; ...
               1 1 1 1 1 1; ...
               3 3 3 3 3 3; ...
               3 3 3 3 3 3; ...
               3 3 3 3 3 3];
    L0      = L0(I(1:nel,1),:);
else
    switch SETTINGS.refinement
        case 'regular'
            L0 = SETTINGS.h0*ones(size(barmid,1),1); % desired length for each bar (spring)
        case 'interface'
            L0 = bar_L0_interface(barmid,SETTINGS.h0,INTERFACE);
        case 'guide_mesh'
            L0 = bar_L0_guide(barmid,GUIDE_MESH); % desired length for each bar (spring)
            
%             % Desired length following Anderson et al. (2005) --> slower and the
%             % results are worse than using bar_L0_guide for barmids.
%             GCOORD    = GCOORD';
%             VCOORD    = GCOORD(:,EL2NOD'); % element node vertex coordinates
%             L0_vertex = reshape(bar_L0_guide(VCOORD',GUIDE_MESH),4,[]); % desired length for each vertex
%             V         = 2./(3.*sum(sqrt(2)./L0_vertex.^3)); % optimal tetrahedral volume
%             L0        = zeros(size(bars,1),1);
%             for i=1:size(bars,1)
%                 EL2bar = find(sum(ismember(EL2NOD,bars(i,:)),2) == 2)';
%                 L0(i)  = (6*sqrt(2))^(1/3).*(sum(V(EL2bar))/size(EL2bar,2)).^(1/3); % desired length for each edge
%             end
    end
    L0      = reshape(L0,nel,6); % desired length of the bars for each element
end
L0_mean     = sum(L0,2)/6;      % average of desired length of the bars for each element
r0_mean     = L0_mean/sqrt(24); % desired inradius for each element (for a regular tetrahedron of side L0_mean)
c           = 1; % constant to relate forces with inradius change (because is not just Hooke's law)
k           = 2; % stiffness (spring constant)
F_tan       = -c*k*(r_I-r0_mean); % ballon force in radial direction applied to either expand or 
                                  % compress the element at each tangent node for each element 
                                  % F>0 implies forces acting outwards the element
                                  % F<0 implies forces acting inwards the element
F_tan       = reshape(repmat(F_tan,1,4)',4*nel,1); % repeat the vector F in order to apply the same forces 
                                                   % to each tangent point for each element (4 tangent points)

%==========================================================================
% COMPUTE TANGENT POINTS (x_tan,y_tan,z_tan) IN WHICH THOSE FORCES ARE
% ACTING
%==========================================================================
x_I         = reshape(repmat(x_I,1,4)',4*nel,1); % every 4 rows is for 1 element (4 faces each element)
y_I         = reshape(repmat(y_I,1,4)',4*nel,1); % every 4 rows is for 1 element (4 faces each element)
z_I         = reshape(repmat(z_I,1,4)',4*nel,1); % every 4 rows is for 1 element (4 faces each element)
r_I         = reshape(repmat(r_I,1,4)',4*nel,1); % every 4 rows is for 1 element (4 faces each element)
n           = zeros(4*nel,3);
for i = 1:nel
    n(4*i-3:4*i,:) = [n1(i,:); n2(i,:); n3(i,:); n4(i,:)]; % every 4 rows is for 1 element (4 faces each element)
end
x_tan        = x_I + r_I.*n(:,1); % x-coordinate for the tangent points, every 4 rows is for 1 element (4 faces each element)
y_tan        = y_I + r_I.*n(:,2); % y-coordinate for the tangent points, every 4 rows is for 1 element (4 faces each element)
z_tan        = z_I + r_I.*n(:,3); % z-coordinate for the tangent points, every 4 rows is for 1 element (4 faces each element)
% data for plotting
PLOT.x_I     = x_I;
PLOT.y_I     = y_I;
PLOT.z_I     = z_I;
PLOT.r_I     = r_I;
PLOT.x_tan   = x_tan;
PLOT.y_tan   = y_tan;
PLOT.z_tan   = z_tan;

%==========================================================================
% ADD NORMAL FORCES TO TANGENT POINTS
%==========================================================================
Fx           = F_tan.*n(:,1);   % x-component of the forces, every 4 rows is for 1 element (4 faces each element)
Fy           = F_tan.*n(:,2);   % y-component of the forces, every 4 rows is for 1 element (4 faces each element)
Fz           = F_tan.*n(:,3);   % z-component of the forces, every 4 rows is for 1 element (4 faces each element)
% data for plotting
PLOT.Fx_end  = x_tan + Fx; % x-coordinate of the end of each force vector, every 4 rows is for 1 element (4 faces each element)
PLOT.Fy_end  = y_tan + Fy; % y-coordinate of the end of each force vector, every 4 rows is for 1 element (4 faces each element)
PLOT.Fz_end  = z_tan + Fz; % z-coordinate of the end of each force vector, every 4 rows is for 1 element (4 faces each element)

%==========================================================================
% SPLIT NORMAL FORCE IN EACH FACE WEIGHING IT BY THE RELATIVE DISTANCE TO 
% THE NODES THAT DEFINE EACH FACE (USING SHAPE FUNCTIONS)
%==========================================================================
N            = shape_function(GCOORD,EL2NOD_sorted,x_tan,y_tan,z_tan);
F_node       = repmat(F_tan,1,4).*N; % Forces (modulus) acting at each node of each face
% every 4 rows is for 1 element (4 faces each element)
F1x          = F_node(:,1).*n(:,1);  % x-component of the forces acting in the 1st node for each element (contribution of 4 faces per element)
F1y          = F_node(:,1).*n(:,2);  % y-component of the forces acting in the 1st node for each element (contribution of 4 faces per element)
F1z          = F_node(:,1).*n(:,3);  % z-component of the forces acting in the 1st node for each element (contribution of 4 faces per element)
F2x          = F_node(:,2).*n(:,1);  % x-component of the forces acting in the 2nd node for each element (contribution of 4 faces per element)
F2y          = F_node(:,2).*n(:,2);  % y-component of the forces acting in the 2nd node for each element (contribution of 4 faces per element)
F2z          = F_node(:,2).*n(:,3);  % z-component of the forces acting in the 2nd node for each element (contribution of 4 faces per element)
F3x          = F_node(:,3).*n(:,1);  % x-component of the forces acting in the 3rd node for each element (contribution of 4 faces per element)
F3y          = F_node(:,3).*n(:,2);  % y-component of the forces acting in the 3rd node for each element (contribution of 4 faces per element)
F3z          = F_node(:,3).*n(:,3);  % z-component of the forces acting in the 3rd node for each element (contribution of 4 faces per element)
F4x          = F_node(:,4).*n(:,1);  % x-component of the forces acting in the 4th node for each element (contribution of 4 faces per element)
F4y          = F_node(:,4).*n(:,2);  % y-component of the forces acting in the 4th node for each element (contribution of 4 faces per element)
F4z          = F_node(:,4).*n(:,3);  % z-component of the forces acting in the 4th node for each element (contribution of 4 faces per element)
face2node    = zeros(4*nel,4);
for i = 1:nel
    face2node(4*i-3:4*i,:) = repmat(EL2NOD_sorted(i,:),4,1);
end
% data for plotting (every 4 rows is for 1 element (4 faces each element))
PLOT.F1x_end = GCOORD(face2node(:,1),1) + F1x; % x-coord of the end of the force vectors acting in the 1st node for each element (contribution of 4 faces per element)
PLOT.F1y_end = GCOORD(face2node(:,1),2) + F1y; % y-coord of the end of the force vectors acting in the 1st node for each element (contribution of 4 faces per element)
PLOT.F1z_end = GCOORD(face2node(:,1),3) + F1z; % z-coord of the end of the force vectors acting in the 1st node for each element (contribution of 4 faces per element)
PLOT.F2x_end = GCOORD(face2node(:,2),1) + F2x; % x-coord of the end of the force vectors acting in the 2nd node for each element (contribution of 4 faces per element)
PLOT.F2y_end = GCOORD(face2node(:,2),2) + F2y; % y-coord of the end of the force vectors acting in the 2nd node for each element (contribution of 4 faces per element)
PLOT.F2z_end = GCOORD(face2node(:,2),3) + F2z; % z-coord of the end of the force vectors acting in the 2nd node for each element (contribution of 4 faces per element)
PLOT.F3x_end = GCOORD(face2node(:,3),1) + F3x; % x-coord of the end of the force vectors acting in the 3rd node for each element (contribution of 4 faces per element)
PLOT.F3y_end = GCOORD(face2node(:,3),2) + F3y; % y-coord of the end of the force vectors acting in the 3rd node for each element (contribution of 4 faces per element)
PLOT.F3z_end = GCOORD(face2node(:,3),3) + F3z; % z-coord of the end of the force vectors acting in the 3rd node for each element (contribution of 4 faces per element)
PLOT.F4x_end = GCOORD(face2node(:,4),1) + F4x; % x-coord of the end of the force vectors acting in the 4th node for each element (contribution of 4 faces per element)
PLOT.F4y_end = GCOORD(face2node(:,4),2) + F4y; % y-coord of the end of the force vectors acting in the 4th node for each element (contribution of 4 faces per element)
PLOT.F4z_end = GCOORD(face2node(:,4),3) + F4z; % z-coord of the end of the force vectors acting in the 4th node for each element (contribution of 4 faces per element)

%==========================================================================
% COMPUTE THE RESULTANT FORCE AT EACH NODE ADDING ALL THE FORCE
% CONTRIBUTIONS (F1x,F1y,F1z,F2x,F2y,F2z,F3x,F3y,F3z,F4x,F4y,F4z)
%==========================================================================
F1234x       = [F1x F2x F3x F4x];
F1234y       = [F1y F2y F3y F4y];
F1234z       = [F1z F2z F3z F4z];
F_balloon    = zeros(3*length(GCOORD),1); % pre-allocating memory
for i = 1:length(GCOORD)
    F_balloon(3*i-2) = sum(F1234x(face2node==i));
    F_balloon(3*i-1) = sum(F1234y(face2node==i));
    F_balloon(3*i  ) = sum(F1234z(face2node==i));
end
F_balloon_x  = F_balloon(1:3:end);
F_balloon_y  = F_balloon(2:3:end);
F_balloon_z  = F_balloon(3:3:end);
% data for plotting
PLOT.F_balloon_x_end = GCOORD(:,1) + F_balloon_x;
PLOT.F_balloon_y_end = GCOORD(:,2) + F_balloon_y;
PLOT.F_balloon_z_end = GCOORD(:,3) + F_balloon_z;

if SETTINGS.show_figs
    figure(203)
    clf
    plot_balloon_forces(GCOORD,EL2NOD_sorted,face2node,PLOT)
end
end % END OF FUNCTION balloon_forces

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function N = shape_function(GCOORD,EL2NOD_sorted,x_tan,y_tan,z_tan)
% Usage: N = shape_function(GCOORD,EL2NOD_sorted,x_tan,y_tan,z_tan)
%
% Purpose:
%   Compute shape functions
%
% Input:
%   GCOORD        : [matrix] : coorfinates of the nodes of the mesh
%   EL2NOD_sorted : [matrix] : conectivity matrix for those elements which
%                              q < q_balloon 
%   x_tan         : [vector] : x-coordinate for the tangent points
%   y_tan         : [vector] : y-coordinate for the tangent points
%   z_tan         : [vector] : z-coordinate for the tangent points
%
% Output:
%   N             : [matrix] : shape functions
%
% JMT Jun 2016

nel = size(EL2NOD_sorted,1); % number of tetrahedrons
x   = reshape(GCOORD(EL2NOD_sorted,1),nel,4)';
y   = reshape(GCOORD(EL2NOD_sorted,2),nel,4)';
z   = reshape(GCOORD(EL2NOD_sorted,3),nel,4)';
xp  = reshape(x_tan,4,nel);
yp  = reshape(y_tan,4,nel);
zp  = reshape(z_tan,4,nel);
r   = zeros(4,nel);
s   = zeros(4,nel);
t   = zeros(4,nel);
for i = 1:4 % compute r,s and t for the 4 tangent points of each element (1 tangent point for each face of the element)
    r(i,:) =   (  (x(2,:)-xp(i,:)).*(y(3,:)-yp(i,:)).*(z(4,:)-zp(i,:)) +      ...
                  (x(3,:)-xp(i,:)).*(y(4,:)-yp(i,:)).*(z(2,:)-zp(i,:)) +      ...
                  (x(4,:)-xp(i,:)).*(y(2,:)-yp(i,:)).*(z(3,:)-zp(i,:)) -      ...
                  (z(2,:)-zp(i,:)).*(y(3,:)-yp(i,:)).*(x(4,:)-xp(i,:)) -      ...
                  (z(3,:)-zp(i,:)).*(y(4,:)-yp(i,:)).*(x(2,:)-xp(i,:)) -      ...
                  (z(4,:)-zp(i,:)).*(y(2,:)-yp(i,:)).*(x(3,:)-xp(i,:))   ) ./ ...
               (  (x(2,:)-x(1,:)) .*(y(3,:)-y(1,:)) .*(z(4,:)-z(1,:))  +      ...
                  (x(3,:)-x(1,:)) .*(y(4,:)-y(1,:)) .*(z(2,:)-z(1,:))  +      ...
                  (x(4,:)-x(1,:)) .*(y(2,:)-y(1,:)) .*(z(3,:)-z(1,:))  -      ...
                  (z(2,:)-z(1,:)) .*(y(3,:)-y(1,:)) .*(x(4,:)-x(1,:))  -      ...
                  (z(3,:)-z(1,:)) .*(y(4,:)-y(1,:)) .*(x(2,:)-x(1,:))  -      ...
                  (z(4,:)-z(1,:)) .*(y(2,:)-y(1,:)) .*(x(3,:)-x(1,:))    );
              
    s(i,:) =   (  (x(1,:)-xp(i,:)).*(y(3,:)-yp(i,:)).*(z(4,:)-zp(i,:)) +      ...
                  (x(3,:)-xp(i,:)).*(y(4,:)-yp(i,:)).*(z(1,:)-zp(i,:)) +      ...
                  (x(4,:)-xp(i,:)).*(y(1,:)-yp(i,:)).*(z(3,:)-zp(i,:)) -      ...
                  (z(1,:)-zp(i,:)).*(y(3,:)-yp(i,:)).*(x(4,:)-xp(i,:)) -      ...
                  (z(3,:)-zp(i,:)).*(y(4,:)-yp(i,:)).*(x(1,:)-xp(i,:)) -      ...
                  (z(4,:)-zp(i,:)).*(y(1,:)-yp(i,:)).*(x(3,:)-xp(i,:))   ) ./ ...
               (  (x(1,:)-x(2,:)) .*(y(3,:)-y(2,:)) .*(z(4,:)-z(2,:))  +      ...
                  (x(3,:)-x(2,:)) .*(y(4,:)-y(2,:)) .*(z(1,:)-z(2,:))  +      ...
                  (x(4,:)-x(2,:)) .*(y(1,:)-y(2,:)) .*(z(3,:)-z(2,:))  -      ...
                  (z(1,:)-z(2,:)) .*(y(3,:)-y(2,:)) .*(x(4,:)-x(2,:))  -      ...
                  (z(3,:)-z(2,:)) .*(y(4,:)-y(2,:)) .*(x(1,:)-x(2,:))  -      ...
                  (z(4,:)-z(2,:)) .*(y(1,:)-y(2,:)) .*(x(3,:)-x(2,:))  );
             
    t(i,:) =   (  (x(1,:)-xp(i,:)).*(y(2,:)-yp(i,:)).*(z(4,:)-zp(i,:)) +      ...
                  (x(2,:)-xp(i,:)).*(y(4,:)-yp(i,:)).*(z(1,:)-zp(i,:)) +      ...
                  (x(4,:)-xp(i,:)).*(y(1,:)-yp(i,:)).*(z(2,:)-zp(i,:)) -      ...
                  (z(1,:)-zp(i,:)).*(y(2,:)-yp(i,:)).*(x(4,:)-xp(i,:)) -      ...
                  (z(2,:)-zp(i,:)).*(y(4,:)-yp(i,:)).*(x(1,:)-xp(i,:)) -      ...
                  (z(4,:)-zp(i,:)).*(y(1,:)-yp(i,:)).*(x(2,:)-xp(i,:))   ) ./ ...
               (  (x(1,:)-x(3,:)) .*(y(2,:)-y(3,:)) .*(z(4,:)-z(3,:))  +      ...
                  (x(2,:)-x(3,:)) .*(y(4,:)-y(3,:)) .*(z(1,:)-z(3,:))  +      ...
                  (x(4,:)-x(3,:)) .*(y(1,:)-y(3,:)) .*(z(2,:)-z(3,:))  -      ...
                  (z(1,:)-z(3,:)) .*(y(2,:)-y(3,:)) .*(x(4,:)-x(3,:))  -      ...
                  (z(2,:)-z(3,:)) .*(y(4,:)-y(3,:)) .*(x(1,:)-x(3,:))  -      ...
                  (z(4,:)-z(3,:)) .*(y(1,:)-y(3,:)) .*(x(2,:)-x(3,:))  );
end
w = 1-r-s-t;
N = zeros(4*nel,4);
for i = 1:nel
    % Shape functions. In columns N1, N2, N3 and N4. Every 4 rows is for 1 element (4 faces each element)
    N(4*i-3:4*i,:) = [r(:,i) s(:,i) t(:,i) w(:,i)];
end
N(N>1)     = 1; % To avoid strange interpolation values due to round off errors
N(N<1e-15) = 0; % To avoid strange interpolation values due to round off errors

end % END OF SUBFUNCTION shape_function

function plot_balloon_forces(GCOORD,EL2NOD_sorted,face2node,PLOT)
% Usage: plot_balloon_forces(GCOORD,EL2NOD_sorted,face2node,PLOT)
%
% Purpose:
%   plot balloon forces
%
% Input:
%   GCOORD        : [matrix]   : coorfinates of the nodes of the mesh
%   EL2NOD_sorted : [matrix]   : conectivity matrix for those elements which
%                                q < q_balloon
%   PLOT          : [structure]: structure containing variables to plot
%
% Output:
%   none (plot figure)
%
% JMT Jun 2016

% Load variables to plot
x_I             = PLOT.x_I;
y_I             = PLOT.y_I;
z_I             = PLOT.z_I;
r_I             = PLOT.r_I;
x_tan           = PLOT.x_tan;
y_tan           = PLOT.y_tan;
z_tan           = PLOT.z_tan;
Fx_end          = PLOT.Fx_end;
Fy_end          = PLOT.Fy_end;
Fz_end          = PLOT.Fz_end;
F1x_end         = PLOT.F1x_end;
F1y_end         = PLOT.F1y_end;
F1z_end         = PLOT.F1z_end;
F2x_end         = PLOT.F2x_end;
F2y_end         = PLOT.F2y_end;
F2z_end         = PLOT.F2z_end;
F3x_end         = PLOT.F3x_end;
F3y_end         = PLOT.F3y_end;
F3z_end         = PLOT.F3z_end;
F4x_end         = PLOT.F4x_end;
F4y_end         = PLOT.F4y_end;
F4z_end         = PLOT.F4z_end;
F_balloon_x_end = PLOT.F_balloon_x_end;
F_balloon_y_end = PLOT.F_balloon_y_end;
F_balloon_z_end = PLOT.F_balloon_z_end;
nel             = size(EL2NOD_sorted,1); % number of tetrahedrons
faceColor = [0.6875 0.8750 0.8984];
tetramesh(EL2NOD_sorted,GCOORD,'FaceColor',faceColor,'FaceAlpha',0.3);
hold on
% text(GCOORD(:,1)+0.1,GCOORD(:,2)+0.1,GCOORD(:,3)+0.1,num2str([1:length(GCOORD)]'),'Color',[0 0 1],'Fontsize',10,'FontWeight','bold')
% scatter3(GCOORD(:,1),GCOORD(:,2),GCOORD(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
scatter3(x_I,y_I,z_I,'rx') % incentre
lightGrey = 0.9*[1 1 1];
for i=1:nel
    [x_sphere,y_sphere,z_sphere] = sphere(50);
    x_sphere = x_sphere*r_I(4*i-3) + x_I(4*i-3);
    y_sphere = y_sphere*r_I(4*i-3) + y_I(4*i-3);
    z_sphere = z_sphere*r_I(4*i-3) + z_I(4*i-3);
    surface(x_sphere,y_sphere,z_sphere,'FaceColor', 'none','EdgeColor',lightGrey)
    hold on
end
% tangent points to the incircle
scatter3(x_tan(1:4:end),y_tan(1:4:end),z_tan(1:4:end),'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]) % face 1 --> green
scatter3(x_tan(2:4:end),y_tan(2:4:end),z_tan(2:4:end),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) % face 2 --> blue
scatter3(x_tan(3:4:end),y_tan(3:4:end),z_tan(3:4:end),'MarkerEdgeColor','k','MarkerFaceColor',[1 1 0]) % face 3 --> yellow
scatter3(x_tan(4:4:end),y_tan(4:4:end),z_tan(4:4:end),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1]) % face 4 --> magenta
% end of normal vector to tangent point at each face
scatter3(Fx_end(1:4:end),Fy_end(1:4:end),Fz_end(1:4:end),'g^') % face 1 --> green
scatter3(Fx_end(2:4:end),Fy_end(2:4:end),Fz_end(2:4:end),'b^') % face 2 --> blue
scatter3(Fx_end(3:4:end),Fy_end(3:4:end),Fz_end(3:4:end),'y^') % face 3 --> yellow
scatter3(Fx_end(4:4:end),Fy_end(4:4:end),Fz_end(4:4:end),'m^') % face 4 --> magenta
% line to link tangent points to the incircle and the end of normal vector (to tangent point at each face)
line([x_tan(1:4:end) Fx_end(1:4:end)]', ...
     [y_tan(1:4:end) Fy_end(1:4:end)]', ...
     [z_tan(1:4:end) Fz_end(1:4:end)]','Color',[0 1 0]) % face 1 --> green
line([x_tan(2:4:end) Fx_end(2:4:end)]', ...
     [y_tan(2:4:end) Fy_end(2:4:end)]', ...
     [z_tan(2:4:end) Fz_end(2:4:end)]','Color',[0 0 1]) % face 2 --> blue
line([x_tan(3:4:end) Fx_end(3:4:end)]', ...
     [y_tan(3:4:end) Fy_end(3:4:end)]', ...
     [z_tan(3:4:end) Fz_end(3:4:end)]','Color',[1 1 0]) % face 3 --> yellow
line([x_tan(4:4:end) Fx_end(4:4:end)]', ...
     [y_tan(4:4:end) Fy_end(4:4:end)]', ...
     [z_tan(4:4:end) Fz_end(4:4:end)]','Color',[1 0 1]) % face 4 --> magenta

scatter3(F1x_end(1:4:end),F1y_end(1:4:end),F1z_end(1:4:end),'gd') % face 1 --> green
scatter3(F2x_end(1:4:end),F2y_end(1:4:end),F2z_end(1:4:end),'gd') % face 1 --> green
scatter3(F3x_end(1:4:end),F3y_end(1:4:end),F3z_end(1:4:end),'gd') % face 1 --> green
scatter3(F4x_end(1:4:end),F4y_end(1:4:end),F4z_end(1:4:end),'gd') % face 1 --> green
line([GCOORD(face2node(1:4:end,1),1) F1x_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,1),2) F1y_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,1),3) F1z_end(1:4:end)]','Color',[0 1 0]) % face 1 --> green
line([GCOORD(face2node(1:4:end,2),1) F2x_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,2),2) F2y_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,2),3) F2z_end(1:4:end)]','Color',[0 1 0]) % face 1 --> green
line([GCOORD(face2node(1:4:end,3),1) F3x_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,3),2) F3y_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,3),3) F3z_end(1:4:end)]','Color',[0 1 0]) % face 1 --> green
line([GCOORD(face2node(1:4:end,4),1) F4x_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,4),2) F4y_end(1:4:end)]', ...
     [GCOORD(face2node(1:4:end,4),3) F4z_end(1:4:end)]','Color',[0 1 0]) % face 1 --> green

scatter3(F1x_end(2:4:end),F1y_end(2:4:end),F1z_end(2:4:end),'bd') % face 2 --> blue
scatter3(F2x_end(2:4:end),F2y_end(2:4:end),F2z_end(2:4:end),'bd') % face 2 --> blue
scatter3(F3x_end(2:4:end),F3y_end(2:4:end),F3z_end(2:4:end),'bd') % face 2 --> blue
scatter3(F4x_end(2:4:end),F4y_end(2:4:end),F4z_end(2:4:end),'bd') % face 2 --> blue
line([GCOORD(face2node(2:4:end,1),1) F1x_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,1),2) F1y_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,1),3) F1z_end(2:4:end)]','Color',[0 0 1]) % face 2 --> blue
line([GCOORD(face2node(2:4:end,2),1) F2x_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,2),2) F2y_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,2),3) F2z_end(2:4:end)]','Color',[0 0 1]) % face 2 --> blue
line([GCOORD(face2node(2:4:end,3),1) F3x_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,3),2) F3y_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,3),3) F3z_end(2:4:end)]','Color',[0 0 1]) % face 2 --> blue
line([GCOORD(face2node(2:4:end,4),1) F4x_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,4),2) F4y_end(2:4:end)]', ...
     [GCOORD(face2node(2:4:end,4),3) F4z_end(2:4:end)]','Color',[0 0 1]) % face 2 --> blue

scatter3(F1x_end(3:4:end),F1y_end(3:4:end),F1z_end(3:4:end),'yd') % face 3 --> yellow
scatter3(F2x_end(3:4:end),F2y_end(3:4:end),F2z_end(3:4:end),'yd') % face 3 --> yellow
scatter3(F3x_end(3:4:end),F3y_end(3:4:end),F3z_end(3:4:end),'yd') % face 3 --> yellow
scatter3(F4x_end(3:4:end),F4y_end(3:4:end),F4z_end(3:4:end),'yd') % face 3 --> yellow
line([GCOORD(face2node(3:4:end,1),1) F1x_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,1),2) F1y_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,1),3) F1z_end(3:4:end)]','Color',[1 1 0]) % face 3 --> yellow
line([GCOORD(face2node(3:4:end,2),1) F2x_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,2),2) F2y_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,2),3) F2z_end(3:4:end)]','Color',[1 1 0]) % face 3 --> yellow
line([GCOORD(face2node(3:4:end,3),1) F3x_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,3),2) F3y_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,3),3) F3z_end(3:4:end)]','Color',[1 1 0]) % face 3 --> yellow
line([GCOORD(face2node(3:4:end,4),1) F4x_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,4),2) F4y_end(3:4:end)]', ...
     [GCOORD(face2node(3:4:end,4),3) F4z_end(3:4:end)]','Color',[1 1 0]) % face 3 --> yellow

scatter3(F1x_end(4:4:end),F1y_end(4:4:end),F1z_end(4:4:end),'md') % face 4 --> magenta
scatter3(F2x_end(4:4:end),F2y_end(4:4:end),F2z_end(4:4:end),'md') % face 4 --> magenta
scatter3(F3x_end(4:4:end),F3y_end(4:4:end),F3z_end(4:4:end),'md') % face 4 --> magenta
scatter3(F4x_end(4:4:end),F4y_end(4:4:end),F4z_end(4:4:end),'md') % face 4 --> magenta
line([GCOORD(face2node(4:4:end,1),1) F1x_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,1),2) F1y_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,1),3) F1z_end(4:4:end)]','Color',[1 0 1]) % face 4 --> magenta
line([GCOORD(face2node(4:4:end,2),1) F2x_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,2),2) F2y_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,2),3) F2z_end(4:4:end)]','Color',[1 0 1]) % face 4 --> magenta
line([GCOORD(face2node(4:4:end,3),1) F3x_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,3),2) F3y_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,3),3) F3z_end(4:4:end)]','Color',[1 0 1]) % face 4 --> magenta
line([GCOORD(face2node(4:4:end,4),1) F4x_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,4),2) F4y_end(4:4:end)]', ...
     [GCOORD(face2node(4:4:end,4),3) F4z_end(4:4:end)]','Color',[1 0 1]) % face 4 --> magenta

scatter3(F_balloon_x_end(EL2NOD_sorted(1,:)),F_balloon_y_end(EL2NOD_sorted(1,:)),F_balloon_z_end(EL2NOD_sorted(1,:)),'r^')
line([GCOORD(:,1) F_balloon_x_end]',[GCOORD(:,2) F_balloon_y_end]',[GCOORD(:,3) F_balloon_z_end]','Color',[1 0 0])
% axis([-6371 6371 -6371 6371 -6371 6371])
axis equal
view(142.5,30)
grid on
xlabel('X');ylabel('Y');zlabel('Z')

end % END OF SUBFUNCTION plot_balloon_forces