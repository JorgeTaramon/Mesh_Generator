function F_balloon = balloon_forces(MESH,EL2NOD_sorted,GUIDE_MESH,SETTINGS)
% Usage: F_balloon = balloon_forces(MESH,EL2NOD_sorted,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   This routine is applied in those elements which quality factor is
%   smaller than the quality tolerance (q < q_balloon).
%   Compute the needed force at each node (for each 'bad' element) taking
%   into account the difference between the incircle for the desired
%   triangle and the incircle for actual triangle. This force will make
%   better elements.
%
% Input:
%   MESH          : [structure]     : structure containing the mesh
%   EL2NOD_sorted : [matrix]        : conectivity matrix for those elements 
%                                     which q < q_balloon 
%   GUIDE_MESH    : [structure]     : structure containing guide mesh     
%   SETTINGS      : [structure]     : structure containing mesh settings
%
% Output:
%   F_balloon     : [column vector] : forces acting at each node of the
%                                     'bad' elements
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if nargin == 0
    MESH.GCOORD        = [0 0; 1 1; 3 0; 3 2; 0 3; 3 3];
    MESH.EL2NOD        = delaunay(MESH.GCOORD);
    q_balloon          = 0.6;
    q                  = mesh_quality(MESH.GCOORD,MESH.EL2NOD);
    [q_sorted,I]       = sort(q);                                % sort q from worst q to best q
    EL2NOD_sorted      = MESH.EL2NOD(I(q_sorted < q_balloon),:); % select those elements which q < q_balloon
    SETTINGS.show_figs = 1;
end

%==========================================================================
% COMPUTE THE INCENTER I(x_I,y_I), THE INRADIUS (ACTUAL,r_I, AND
% DESIRED,r0_mean) AND THE FORCES ACTING AT EACH TANGENT POINT OF EACH
% ELEMENT
%==========================================================================
%                 3
%                 .
%                / \                           [                      ]
%               /   \                          |  a*x1 + b*x2 + c*x3  |
%              /     \                   I_x = | -------------------- |
%             /       \                        |      a + b + c       |
%            /         \                       [                      ]
%         b /           \ a
%          /             \                     [                      ]
%         /   (I_x,I_y)   \                    |  a*y1 + b*y2 + c*y3  |
%        /        *        \             I_y = | -------------------- |
%       /                   \                  |      a + b + c       |
%      /                     \                 [                      ]
%     .-----------------------.
%    1            c            2

GCOORD = MESH.GCOORD;
EL2NOD = MESH.EL2NOD;
nel    = size(EL2NOD_sorted,1);
bars   = [EL2NOD_sorted(:,[1,2]);...
          EL2NOD_sorted(:,[2,3]);...
          EL2NOD_sorted(:,[3,1])];
barvec = GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:);
barmid = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2;
L      = sqrt(sum(barvec.^2,2)); % bar length
L      = reshape(L,size(EL2NOD_sorted,1),3); % bar length for each element
I_x    = ((L(:,1).*GCOORD(EL2NOD_sorted(:,3),1))+...
          (L(:,2).*GCOORD(EL2NOD_sorted(:,1),1))+...
          (L(:,3).*GCOORD(EL2NOD_sorted(:,2),1)))./(sum(L,2));
         % x-coordinate of the incenter for each element
I_y    = ((L(:,1).*GCOORD(EL2NOD_sorted(:,3),2))+...
          (L(:,2).*GCOORD(EL2NOD_sorted(:,1),2))+...
          (L(:,3).*GCOORD(EL2NOD_sorted(:,2),2)))./(sum(L,2));
         % y-coordinate of the incenter for each element
r_I    = sqrt(((L(:,2)+L(:,3)-L(:,1)).*...
               (L(:,3)+L(:,1)-L(:,2)).*...
               (L(:,1)+L(:,2)-L(:,3)))./...
               (L(:,1)+L(:,2)+L(:,3)))/2; % inradius for each element
if nargin == 0
    L0 = [3;1;3;3;1;3;3;1;3];
else
    if strcmp(SETTINGS.guide_mesh,'no')
        L0 = h0*ones(size(barmid,1),1); % desired length for each bar
    else
        L0 = bar_L0_guide(barmid,GUIDE_MESH,SETTINGS); % desired length for each bar
    end
end
L0      = reshape(L0,nel,3); % desired length of the bars for each element
L0_mean = sum(L0,2)/3;         % average of desired length of the bars for each element
r0_mean = (sqrt(3)/6)*L0_mean; % desired inradius for each element
c       = 1; % constant to relate forces with inradius change (because is not just Hooke's law)
k       = 1; % stiffness (spring constant)
F       =  -c*k*(r_I-r0_mean); % balloon force in radial direction applied to either expand or 
                               % compress the element at each tangent node for each element
                               % F > 0 implies forces acting outwards the element
                               % F < 0 implies forces acting inwards the element
F       = repmat(F,3,1); % repeat the vector F in order to have the same elements as in barvec

%==========================================================================
% COMPUTE TANGENT POINTS (x_tan,y_tan) IN WHICH THOSE FORCES ARE ACTING
%==========================================================================
% Compute the straight-line equations that are perpendicular to each bar and pass through the incenters
barvec_normal = [barvec(:,2) -barvec(:,1)]; % perpendicular vector to each bar (outside direction)
barvec_normal = repmat(sign(F),1,2).*barvec_normal; % give orientation to the perpendicular vector depending on the sign of F
                                                    % if F > 0 the perpendicular vector will point outwards the element
                                                    % if F < 0 the perpendicular vector will point inwards the element
I_x   = repmat(I_x,3,1); % repeat the vector I_x in order to have the same elements as in barvec
I_y   = repmat(I_y,3,1); % repeat the vector I_y in order to have the same elements as in barvec
alpha = atan2d(barvec(:,2),barvec(:,1)); % angle that defines the vector for each bar

% Find the intersection point between the bars and the perpendicular lines
% to bars that contain the incenters. The intersection of two lines L1 and
% L2 in 2-D with L1 containting the points (x1,y1) and (x2,y2), and L2 
% containing the points (x3,y3)and (x4,y4) is given by:
%
%         | |x1  y1|             |               | |x1  y1|             |
%         | |      |     x1 - x2 |               | |      |     y1 - y2 |
%         | |x2  y2|             |               | |x2  y2|             |
%         |                      |               |                      |
%         |                      |               |                      |
%         | |x3  y3|             |               | |x3  y3|             |
%         | |      |     x3 - x4 |               | |      |     y3 - y4 |
%         | |x4  y4|             |               | |x4  y4|             |
%    x = --------------------------         y = -------------------------- 
%         | x1 - x2      y1 - y2 |               | x1 - x2      y1 - y2 |
%         |                      |               |                      |
%         | x3 - x4      y3 - y4 |               | x3 - x4      y3 - y4 |

% We have to compute the 4th point since we only have 3: the ends of the
% bars ((x1,y1) and (x2,y2)) and the incenter (x3,y3). For this, we carry
% the distance L from the incenter along the perpendicular line to bar
x1    = GCOORD(bars(:,1),1);
y1    = GCOORD(bars(:,1),2);
x2    = GCOORD(bars(:,2),1);
y2    = GCOORD(bars(:,2),2);
x3    = I_x;
y3    = I_y;
L     = reshape(L,length(F),1); % bar lengths
x4    = x3 + L.*cosd(alpha-90); % x-coordinate for the 4th point
y4    = y3 + L.*sind(alpha-90); % y-coordinate for the 4th point
det12 = x1.*y2 - y1.*x2;
det34 = x3.*y4 - y3.*x4;
x1_x2 = (x1-x2);
x3_x4 = (x3-x4);
y1_y2 = (y1-y2);
y3_y4 = (y3-y4);
x_tan = (det12.*x3_x4 - x1_x2.*det34)./(x1_x2.*y3_y4 - y1_y2.*x3_x4);
y_tan = (det12.*y3_y4 - y1_y2.*det34)./(x1_x2.*y3_y4 - y1_y2.*x3_x4);

%==========================================================================
% ADD NORMAL FORCES TO TANGENT POINTS
%==========================================================================
alpha_normal = atan2d(barvec_normal(:,2),barvec_normal(:,1)); 
% angle that defines the perpendicular vectors(forces) to the tangent points
% data for plotting
Fx           = abs(F).*cosd(alpha_normal); % x-component of the normal forces
Fy           = abs(F).*sind(alpha_normal); % y-component of the normal forces
Fx_end       = x_tan + Fx; % x-coordinate of the end of each normal force vector
Fy_end       = y_tan + Fy; % y-coordinate of the end of each normal force vector

%==========================================================================
% SPLIT NORMAL FORCE IN EACH BAR WEIGHING IT BY THE RELATIVE DISTANCE TO 
% THE ENDS OF THE BAR
%==========================================================================
d      = sqrt((GCOORD(bars(:,1),1)-x_tan).^2+...
              (GCOORD(bars(:,1),2)-y_tan).^2); % distance from tangent point to the first nodes in bars
F1     = abs(F).*((L-d)./L); % Force (modulus) acting in the 1st bar nodes. 
                             % The force is weighed by the distance from
                             % the tangent point to the 1st bar node.
F2     = abs(F).*(d./L);     % Force (modulus) acting in the 2nd bar nodes.
                             % The force is weighed by the distance from
                             % the tangent point to the 2nd bar node 
F1x   = abs(F1).*cosd(alpha_normal); % x-component of the force that acts in the 1st bar nodes
F1y   = abs(F1).*sind(alpha_normal); % y-component of the force that acts in the 1st bar nodes
F2x   = abs(F2).*cosd(alpha_normal); % x-component of the force that acts in the 2nd bar nodes
F2y   = abs(F2).*sind(alpha_normal); % y-component of the force that acts in the 2nd bar nodes
% data for plotting
F1x_end = GCOORD(bars(:,1),1) + F1x; % x-coordinate of the end of the force vectors acting in the 1st bar nodes
F1y_end = GCOORD(bars(:,1),2) + F1y; % y-coordinate of the end of the force vectors acting in the 1st bar nodes
F2x_end = GCOORD(bars(:,2),1) + F2x; % x-coordinate of the end of the force vectors acting in the 2nd bar nodes
F2y_end = GCOORD(bars(:,2),2) + F2y; % y-coordinate of the end of the force vectors acting in the 2nd bar nodes

%==========================================================================
% COMPUTE THE RESULTANT FORCE AT EACH NODE ADDING ALL THE FORCE
% CONTRIBUTIONS (F1x,F1y,F2x,F2y)
%==========================================================================
F12x      = [F1x F2x];
F12y      = [F1y F2y];
F_balloon = zeros(2*length(GCOORD),1); % pre-allocating memory
for i=1:length(GCOORD)
    F_balloon(2*i-1) = sum(F12x(bars==i));
    F_balloon(2*i)   = sum(F12y(bars==i));
end
% data for plotting
F_balloon_x     = F_balloon(1:2:end);
F_balloon_y     = F_balloon(2:2:end);
F_balloon_x_end = GCOORD(F_balloon_x ~= 0,1) + F_balloon_x(F_balloon_x ~= 0); % x-coordinate of the end of the resultant force vector
F_balloon_y_end = GCOORD(F_balloon_y ~= 0,2) + F_balloon_y(F_balloon_y ~= 0); % y-coordinate of the end of the resultant force vector

if SETTINGS.show_figs
    figure(203)
    clf
    trimesh(EL2NOD,GCOORD(:,1),GCOORD(:,2),'Color',[0.5 0.5 0.5])
    hold on
%     text(GCOORD(:,1)+0.1,GCOORD(:,2)+0.1,num2str([1:length(GCOORD)]'),'Color',[0 0 1],'Fontsize',10,'FontWeight','bold')
    plot(x_tan,y_tan,'go')     % tangent points to the incircle
    plot(I_x,I_y,'rx')         % incentre
    plot(Fx_end,Fy_end,'g^')   % end of normal vector to each edge
    plot(F1x_end,F1y_end,'r^') % end of normal vector to each edge at one of the ends (weighted by the distance)
    plot(F2x_end,F2y_end,'r^') % end of normal vector to each edge at the other end (weighted by the distance)
    plot(F_balloon_x_end,F_balloon_y_end,'b^') % end of resultant force vector
    for i=1:size(EL2NOD_sorted,1)
        plot_circle(I_x(i),I_y(i),r_I(i))
    end
    line([x_tan Fx_end]',[y_tan Fy_end]','Color',[0 1 0])
    line([GCOORD(bars(:,1),1) F1x_end]',[GCOORD(bars(:,1),2) F1y_end]','Color',[1 0 0])
    line([GCOORD(bars(:,2),1) F2x_end]',[GCOORD(bars(:,2),2) F2y_end]','Color',[1 0 0])
    line([GCOORD(F_balloon_x ~= 0,1) F_balloon_x_end]',[GCOORD(F_balloon_y ~= 0,2) F_balloon_y_end]','Color',[0 0 1])
    axis equal
end
end % END OF FUNCTION balloon_forces

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function plot_circle(x,y,r)
% Usage: plot_circle(x,y,r)
%
% Purpose:
%   Plot a circle. The angle step is 0.01, bigger values will draw the
%   circle faster but not very smooth.
%
% Input:
%   x : [scalar] : x-coordinate of the circle center
%   y : [scalar] : y-coordinate of the circle center
%   r : [scalar] : radius of the circle 
%
% Output:
%   none (plot figure)
%
% JMT May 2016

ang = 0:0.1:2*pi;
xp  = r*cos(ang);
yp  = r*sin(ang);
plot(x+xp,y+yp);

end % END OF SUBFUNCTION plot_circle