function [L0,L,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
% Usage: [L0,L,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,INTERFACE,SETTINGS)
%
% Purpose:
%   Compute the desired and actual bar lengths of the mesh, the bars, the
%   number of bars and the bar vectors.
%
% Input:
%   MESH       : [structure]     : structure containing the mesh
%   GUIDE_MESH : [structure]     : structure containing guide mesh
%   INTERFACE  : [structure]     : structure containing interface settings
%   SETTINGS   : [structure]     : structure containing mesh settings
%
% Output:
%   L0         : [column vector] : desired bar length
%   L          : [column vector] : actual bar length
%   bars       : [matrix]        : bars
%   nbar       : [scalar]        : number of bars
%   barvec     : [matrix]        : bar vectors
%
% JMT May 2016
% JMT Jan 2017: added sph_interface case to refine the mesh
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

GCOORD     = MESH.GCOORD;
nfix       = size(MESH.pfix,1);
EL2NOD     = MESH.EL2NOD;
bars       = [EL2NOD(:,[1,2]); ...
              EL2NOD(:,[1,3]); ...
              EL2NOD(:,[1,4]); ...
              EL2NOD(:,[2,3]); ...
              EL2NOD(:,[2,4]); ...
              EL2NOD(:,[3,4])];
[bars,~,~] = unique(sort(bars,2),'rows');
if nfix > 12
    % we cannot modify bars whose ends are fixed nodes
    bars(sum(ismember(bars,(1:nfix)),2)==2,:) = []; % remove bars formed only by fixed nodes
end
nbar       = size(bars,1);                                  % number of bars
barvec     = GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:);     % list of bar vectors
L          = sqrt(sum(barvec.^2,2));                        % bar lengths
barmid     = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoint of each bar
switch SETTINGS.refinement
    case 'regular'
        L0 = SETTINGS.h0*ones(size(barmid,1),1);
    case 'interface'
        L0 = bar_L0_interface(barmid,SETTINGS.h0,INTERFACE);
    case 'guide_mesh'
        L0 = bar_L0_guide(barmid,GUIDE_MESH);
        
%         % Desired length following Anderson et al. (2005) --> slower and the
%         % results are worse than using bar_L0_guide for barmids.
%         GCOORD    = GCOORD';
%         VCOORD    = GCOORD(:,EL2NOD'); % element node vertex coordinates
%         L0_vertex = reshape(bar_L0_guide(VCOORD',GUIDE_MESH),4,[]); % desired length for each vertex
%         V         = 2./(3.*sum(sqrt(2)./L0_vertex.^3)); % optimal tetrahedral volume
%         L0        = zeros(size(bars,1),1);
%         for i=1:size(bars,1)
%             EL2bar = find(sum(ismember(EL2NOD,bars(i,:)),2) == 2)';
%             L0(i)  = (6*sqrt(2))^(1/3).*(sum(V(EL2bar))/size(EL2bar,2)).^(1/3);
%         end
end

end % END OF FUNCTION bar_length