function [xyz,tri] = sph_shell3(nlevel,a)
% This function generates a set of points on a sphere of radius a (INPUT)
% by recursively splitting an initial dodecahedron of points nlevel (INPUT) times
%
% INPUT a --- radius of spherical shell
%       nlevel -- number of times to split an initial dodecahedron of 12 vertices & 20 tri
%
% OUTPUT xyz --- x,y,z,coords of the points on the sphere
%        tri --- 6-node triangle connectivity matrix for these points
%
% NOTE The initial polar nodes are 89.9 and -89.9 latitude so azimuths are
% definable
% 
%Based on sphgen2.f FORTRAN code
% c...JPM&WJM 7/8/92   mod. 7 April 98 // 26 April 99, 6Sept 99
% c
% c         level   #triangles    #nodes	typ.angular dist (times 111.11 for km)
% c            1           20         42	  31.68
% c            2           80        162	  15.79
% c            3          320        642	   7.96
% c            4         1280       2562	   3.98
% c            5         5120      10242	   1.99
% c            6        20480      40962	    .99
% c            7        81920     163842        .50
% c            8       327680     655362        .25
% c            9      1310720    2621442        .13
%
%
%
if nargin == 0
    nlevel = 1; % if no input => TEST
    a = 1.0;
end

%  SET UP INITIAL ICOSAHEDRON -- tri is connectivity matrix for triangles,
%                                xyz are xyz coordinates of node-points (vectors on sphere)
%xyz = zeros(3,12);
[xyz,tri(:,:),nnodes,ntri] = icosahedron(); % 20 starting triangles (icosahedron)
%
%  ADD MIDSIDE NODES 4-6 for each icosahedron triangle
%
% new vectorized way, removing duplicates with unique_keep_order
sides = reshape(tri([1,2,2,3,3,1],:),2,[])'; % sides of each tri, in order 1,2,3 (trinodes 4,5,6)
[sorted_sides,JJ] = unique_keep_order(sort(sides,2)); % sort and remove duplicate sides
     % returns unique sides of sorted_sides=sort(sides,2) which sorts side-vertex pairs
     % in the order they first appear in the list JJ => sorted_sides(:) = sides(JJ) USED BELOW    
     % Use sort(sides,2) in rhs because you want all side node-pairs to be in same order 
     % to find the pair of matching sides for two neighbor tri's (low node num, high node num)
     % e.g.. side connecting nodes (2,1) is same side as (1,2) - after sort (2,1) would be (1,2)
% calculate location of new side nodes (incl doubles)
xyz_side = xyz(:,sorted_sides(:,1)) + xyz(:,sorted_sides(:,2));
% make them vector length 1
scale    = sum( xyz_side.^2 ); % dot product -- unrolled loop for l^2
xyz_side = xyz_side ./ repmat(sqrt(scale),[3 1]); % scale them all by l
                                % repmat 3 copies of l, 1 for each component
% cat side nodes behind vertex nodes
xyz      = [xyz xyz_side];
    
nsides              = length(sorted_sides);
node_number_side    = [1:nsides]' + nnodes;
sorted_node_number  = node_number_side(JJ); % Use order JJ that they first appear in sort list 
tri(4:6,:)          = reshape(sorted_node_number,3,[]); % fill in 4:6 of each tri connectivity column
    
nnodes   = size(xyz,2);
nvertices = nnodes-nsides; % used in plotting

%
% MAIN REFINEMENT LOOP
%
for ilevel = 2:nlevel
     tri_old = tri;       % start from previous mesh of triangles
     ntri    = 4*ntri;    % each old triangle will be split into 4 new triangles
     tri        = zeros(6,ntri);% zero array for new level of 6 node triangles
     tri(1:3,:) = reshape(tri_old([1 4 6; 4 2 5; 6 5 3; 4 5 6]',:),3,[]);% vertex nodes of 4 new tri's
    % ---------------------------------------------------------------------
    % FORTRAN loop over triangle including all functions called from within "genmid" is
    % replaced by the next lines. JH & JPM 16April2010

    sides = reshape(tri([1,2,2,3,3,1],:),2,[])'; % sides of each triangle, in order 1,2,3 (will be nodes 4,5,6)
    [sorted_sides,JJ] = unique_keep_order(sort(sides,2)); % sort and remove duplicate sides
           % returns unique sides of sorted_sides=sort(sides,2) which sorts side-vertex pairs
           % in the order they first appear in the list JJ => sorted_sides(:) = sides(JJ) USED BELOW    
           % Use sort(sides,2) in rhs because you want all side node-pairs to be in same order 
           % to find the pair of matching sides for two neighbor tri's (low node num, high node num)
           % e.g.. side connecting nodes (2,1) is same side as (1,2) - after sort (2,1) would be (1,2)
    % calculate location of new side nodes (incl doubles)
    xyz_side = xyz(:,sorted_sides(:,1)) + xyz(:,sorted_sides(:,2));
    % make them vector length 1
    scale    = sum( xyz_side.^2 ); % dot product -- unrolled loop for l^2 vectorized to work on all pts
    xyz_side = xyz_side ./ repmat(sqrt(scale),[3 1]); % scale them all by l
                                % repmat 3 copies of l, 1 for each component
    % cat side nodes behind vertex nodes
    xyz      = [xyz xyz_side];
    
    nsides             = length(sorted_sides);
    node_number_side   = [1:nsides]' + nnodes;
    sorted_node_number = node_number_side(JJ);
    tri(4:6,:)         = reshape(sorted_node_number,3,[]); % make columns
    
    nnodes   = size(xyz,2);
    nvertices = nnodes-nsides; % used in plotting

    
    % ---------------------------------------------------------------------
end
%----------------------------
% END OF MAIN REFINEMENT LOOP
%----------------------------

% scale xyz by radius a (right now is a unit sphere)
%
xyz = a .* xyz;

%
% % plot mesh
% figure(2);clf;hold on
% h = trimesh(tri(1:3,:)',xyz(1,:)',xyz(2,:)',xyz(3,:)');
% set(h,'edgecolor','k');
% plot3(xyz(1,1:nvertices)',xyz(2,1:nvertices)',xyz(3,1:nvertices)','o','MarkerFaceColor','b','MarkerSize',7); %plot xyz vertex points
% plot3(xyz(1,(nvertices+1):nnodes)',xyz(2,(nvertices+1):nnodes)',xyz(3,(nvertices+1):nnodes)','ob','MarkerSize',6); % plot of xyz midside points (red)
% if nlevel < 2
%     %vertex nodes -- 14 pt
%     text(xyz(1,1:nvertices)',xyz(2,1:nvertices)',xyz(3,1:nvertices)',...
%          num2str([1:nvertices]'),'FontSize',14,'Color','b');
%     %midside nodes -- 9 pt
%     text(xyz(1,(nvertices+1):nnodes)',xyz(2,(nvertices+1):nnodes)',xyz(3,(nvertices+1):nnodes)',...
%          num2str([(nvertices+1):nnodes]'),'FontSize',9,'Color','r');
% end
% view(-230,40); axis equal;

        %
        % SORT nodes by decreasing lat, increasing long (used in FORTRAN code to reduce matrix bandwidth, 
        % not needed for MATLAB where we can keep original node number for clarity)
        %
        % alatlong = zeros(2,nnodes);
        % distorder = zeros(nnodes,1);
        % for inode = 1:nnodes
        %     [alatlong(1:2,inode)] = vang(xyz(:,inode));    
        % end
        % dist_const = 10*2^nlevel
        % dist_order(1:nnodes,1) = dist_const .* (90 - alatlong(1,1:nnodes)) + alatlong(2,1:nnodes);
        % [wrap_order,IX,IY] = unique(dist_order); % here we want IX,IY -- the permutation indices for the sort, 
        %                                     % instead of wrap_order which is the sorted list of wraparound distances
        %                                     % IX = iold(inew) -- sorted list of old nodes listed in new order
        %                                     % IY = inew(iold) -- new node listed in order corresponding to order of old list 
        %                                   
        % %
        % xyzsort = zeros(size(xyz));                 % define new array for xyz points
        % lmn = zeros(size(tri,1),size(tri,2)); % define new array for lmn connectivity matrix for 6-node triangles
        % %
        % xyzsort(:,:) = xyz(:,IX(:));       % sorted xyz nodes (to be overwritten onto xyz when scaling by radius)
        % lmn = IY(tri(:,:,nlevel));      % sorted connectivity matrix lmn in terms of new distance-sorted node numbers
        % 
        % 
        % %
        % % scale sorted xyz by radius a (right now is a unit sphere)
        % %
        % xyz = a .* xyzsort;
        % % plot3(xyz(1,1:50),xyz(2,1:50),xyz(3,1:50),'*r'); hold on
        % % plot3(xyz(1,51:80),xyz(2,51:80),xyz(3,51:80),'*g'); hold on
        % % plot3(xyz(1,81:nnodes),xyz(2,81:nnodes),xyz(3,81:nnodes),'*b');
        % figure(2);clf;hold on
        % plot3(xyz(1,:)',xyz(2,:)',xyz(3,:)','*'); % simple plot of xyz points
        % text(xyz(1,:)',xyz(2,:)',xyz(3,:)',num2str([1:size(xyz,2)]'))
        % trimesh(lmn(1:3,:)',xyz(1,:)',xyz(2,:)',xyz(3,:)');

end
%===============================================================
% SUBFUNCTIONS USED BY MAIN PROGRAM
%===============================================================
function [xyz,tri,nnodes,ntri] = icosahedron()
%
alatlong = zeros(2,12); % initial lat/long of points of icosahedron
xyz = zeros(3,12);
tri = zeros(6,20); % icosahedron triangles
%c...set up lat, long of icosahedron nodal vertices  (note 'pole' at 89.9)
      alpha = 72.0;
      beta  = 26.56505118;
      alatlong(1,1)  = 89.9999;
      alatlong(2,1)  = 0.0;
      alatlong(1,2)  = beta;
      alatlong(2,2)  = 0.0;
      alatlong(1,3)  = beta;
      alatlong(2,3)  = alpha;
      alatlong(1,4)  = beta;
      alatlong(2,4)  = 2*alpha;
      alatlong(1,5)  = beta;
      alatlong(2,5)  = 3*alpha;
      alatlong(1,6)  = beta;
      alatlong(2,6)  = 4*alpha;
      alatlong(1,7)  =-beta;
      alatlong(2,7)  = 0.5*alpha;
      alatlong(1,8)  =-beta;
      alatlong(2,8)  = 1.5*alpha;
      alatlong(1,9)  =-beta;
      alatlong(2,9)  = 2.5*alpha;
      alatlong(1,10) =-beta;
      alatlong(2,10) = 3.5*alpha;
      alatlong(1,11) =-beta;
      alatlong(2,11) = 4.5*alpha;
      alatlong(1,12) =-89.9999;
      alatlong(2,12) = 0.0;

%c...convert lat-long of nodal vertices to direction cosines
%       for ii = 1:12
%         [xyz(:,ii)] = vset(alatlong(1,ii), alatlong(2,ii) );
%       end
      xyz = vset(alatlong(1,1:12), alatlong(2,1:12) );
%c...set up icosahedron triangles, nodal structure
       
      tri(1,1) = 1;  tri(2,1) = 2;   tri(3,1) = 3;
      tri(1,2) = 1;  tri(2,2) = 3;   tri(3,2) = 4;
      tri(1,3) = 1;  tri(2,3) = 4;   tri(3,3) = 5;
      tri(1,4) = 1;  tri(2,4) = 5;   tri(3,4) = 6;
      tri(1,5) = 1;  tri(2,5) = 6;   tri(3,5) = 2;
      tri(1,6) = 2;  tri(2,6) = 7;   tri(3,6) = 3;
      tri(1,7) = 3;  tri(2,7) = 7;   tri(3,7) = 8;
      tri(1,8) = 3;  tri(2,8) = 8;   tri(3,8) = 4;
      tri(1,9) = 4;  tri(2,9) = 8;   tri(3,9) = 9;
      tri(1,10) = 4; tri(2,10) = 9;  tri(3,10) = 5;
      tri(1,11) = 5; tri(2,11) = 9;  tri(3,11) = 10;
      tri(1,12) = 5; tri(2,12) = 10; tri(3,12) = 6;
      tri(1,13) = 6; tri(2,13) = 10; tri(3,13) = 11;
      tri(1,14) = 6; tri(2,14) = 11; tri(3,14) = 2;
      tri(1,15) = 2; tri(2,15) = 11; tri(3,15) = 7;
      tri(1,16) = 7; tri(2,16) = 12; tri(3,16) = 8;
      tri(1,17) = 8; tri(2,17) = 12; tri(3,17) = 9;
      tri(1,18) = 9; tri(2,18) = 12; tri(3,18) = 10;
      tri(1,19) = 10;tri(2,19) = 12; tri(3,19) = 11;
      tri(1,20) = 11;tri(2,20) = 12; tri(3,20) = 7;

      nnodes = 12;
      ntri = 20;
end
%---------------------------------------------
function [ltri,xyz,nnodes] = genmid(ltri,xyz,nnodes)
%C  GENERATE MID-SIDE NODES FOR A TRIANGLE LMN
%C  (BOTH XYZ NODE VALUES AND ELEMENT CONNECTIVITY)
%C  PROGRAM WILL NOT WORK IF ASKED TO INTERPOLATE TWO EXACTLY OPPOSING POINTS
%C  (NO UNIQUE MIDPOINT IN THIS PATHOLOGICAL CASE WHICH DOES NOT ARISE HERE)
%
%c...node 4 is between nodes 1 and 2
      u(1:3,1) = xyz(:,ltri(1)) + xyz(:,ltri(2));
% convert to unit vector
      u = unitvec(u);
% test if duplicate
      [xyz,MNODE,nnodes] = testdup(u,xyz,nnodes);
      ltri(4)  = MNODE;

%c...node 5 is between nodes 2 and 3
      u(:) = xyz(:,ltri(2)) + xyz(:,ltri(3));
% convert to unit vector
      [u] = unitvec(u);
% test if duplicate
      [xyz,MNODE,nnodes] = testdup(u,xyz,nnodes);
      ltri(5)  = MNODE;

%c...node 6 is between nodes 3 and 1
      u(:) = xyz(:,ltri(3)) + xyz(:,ltri(1));
% convert to unit vector
      [u] = unitvec(u);
% test if duplicate
      [xyz,MNODE,nnodes] = testdup(u,xyz,nnodes);
      ltri(6)  = MNODE;
end
%---------------------------------------------
function [Au,JJu] = unique_keep_order(A)
% function to return unique rows Au of A in the order first encountered
% JJu returns an index to recreate A from Au, i.e. A(:) = Au(JJu)
[tmp,iu,ju] = unique(A,'rows','first'); % find FIRST occurance of unique rows of A
% tmp is an already SORTED version of A, we do NOT want this sort done
[uko1,tmp,JJu] = unique(iu(ju)); % create index uko1 to extract Au from A
                                  % JJ returned to index all of A onto Au
Au = A(uko1,:);
end % END OF FUNCTION unique_keep_order
%---------------------------------------------
function  [xyz,MNODE,nnodes] = testdup(u,xyz,nnodes)
%c...find if node duplicates a previous node
%c...if so return previous node in variable MNODE
%c...if not, add node to list and return MNODE=nnodes + 1
%c...and nnodes = nnodes + 1
%c...where nnodes (input)is number of previous nodes
%c...uses XYZ(3,nnodes) -- list of previous node direction cosines
%c...U is the direction cosines of current point
     criter = .99999; % criterion to define points are same
  for ii = 1:nnodes
        MNODE = ii;
        if(dot(xyz(:,ii),u(:))>criter) 
            return % node does exist, is MNODE
        end
  end
        
%c...fell through loop, 
%c...previous node doesn't exist, add extra node at end of list

   nnodes = nnodes + 1;
   MNODE = nnodes;
   xyz(:,MNODE) = u(:);
end  
%---------------------------------------------
function [newtriangles] = splittri(ltriangle)
%C  SPLIT ONE 6-NODE TRIANGLE INTO FOUR 3-NODE SUB-TRIANGLES (no mid-nodes yet)
newtriangles = zeros(6,4);

newtriangles(1,1) = ltriangle(1);
newtriangles(2,1) = ltriangle(4);
newtriangles(3,1) = ltriangle(6);     

newtriangles(1,2) = ltriangle(4);
newtriangles(2,2) = ltriangle(2);
newtriangles(3,2) = ltriangle(5);    

newtriangles(1,3) = ltriangle(5);
newtriangles(2,3) = ltriangle(6);
newtriangles(3,3) = ltriangle(4);      
 
newtriangles(1,4) = ltriangle(6);
newtriangles(2,4) = ltriangle(5);
newtriangles(3,4) = ltriangle(3);

end
%
%---------------------------------------------
% function [a] = vset(alat,along)
% %c convert  lat,long in degrees into direction cosines  26apr99
%   a = zeros(3,1);
%   deg2rad =  0.0174532925199;
%   a(1) = cos(deg2rad*alat) * cos(deg2rad*along);
%   a(2) = cos(deg2rad*alat) * sin(deg2rad*along);
%   a(3) = sin(deg2rad*alat);
% end
%---------------------------------------------
function [a] = vset(alat,along)
%convert  lat,long in degrees into direction cosines  26apr99
%vectorized  
  a = zeros(3,max(size(alat)));
  deg2rad =  0.0174532925199;
  a(1,:) = cos(deg2rad*alat(:)) .* cos(deg2rad*along(:));
  a(2,:) = cos(deg2rad*alat(:)) .* sin(deg2rad*along(:));
  a(3,:) = sin(deg2rad*alat(:));
end
%---------------------------------------------
function [u] = unitvec(a)
% convert vector into a unit vector
  scale =  dot(a(:),a(:));
  if ( scale <= 0.0 )  
     print  "UNITVEC=0, BLOWUP -- ANSWER NONSENSE"
  else
     u(1:3,1) = a(1:3,1) ./ sqrt(scale);
  end
end
%---------------------------------------------
function [alatlong] = vang(a)
%c        convert from direction cosines into  lat,long   26apr99
%c        (vector 'a' does not have to be of unit length.)
  rad2deg =  1./0.0174532925199;
  aa = sqrt( a(1)*a(1) + a(2)*a(2) );
  alat = 90.*sign( a(3) );
  along = 0.0;
  if (aa > 0.0 ) 
     alat  = atan ( a(3) / aa  ) * rad2deg;
     along = atan2( a(2), a(1) ) * rad2deg;
  end
%c...(uncomment next lines to go from 0-360 degrees instead of -180 to +180)
  if(along < 0.0 ) 
      along = along + 360.;
  end
  if(along > 359.999)
      along = 0.;
  end
  alatlong = [alat;along];
end




	
