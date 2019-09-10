function [L0,L,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,SETTINGS)
% Usage: [L0,L,bars,nbar,barvec] = bar_length(MESH,GUIDE_MESH,SETTINGS)
%
% Purpose:
%   Compute the desired and actual bar lengths of the mesh, the bars, the
%   number of bars and the bar vectors.
%
% Input:
%   MESH       : [structure] : structure containing the mesh
%   GUIDE_MESH : [structure] : structure containing guide mesh
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   L0     : [column vector] : desired bar length
%   L      : [column vector] : actual bar length
%   bars   : [matrix]        : bars
%   nbar   : [scalar]        : number of bars
%   barvec : [matrix]        : bar vectors
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

GCOORD     = MESH.GCOORD;
nfix       = size(MESH.pfix,1);
EL2NOD     = MESH.EL2NOD;
bars       = [EL2NOD(:,[1,2]);EL2NOD(:,[2,3]);EL2NOD(:,[3,1])];
[bars,~,~] = unique(sort(bars,2),'rows');
% if strcmp(SETTINGS.first_guess,'load') && nfix > 12
%     % we cannot modify bars whose ends are fixed nodes
%     bars(sum(ismember(bars,(1:nfix)),2)==2,:) = []; % remove bars formed only by fixed nodes
% end
nbar       = size(bars,1);                                  % number of bars
barvec     = GCOORD(bars(:,2),:) - GCOORD(bars(:,1),:);     % list of bar vectors
L          = sqrt(sum(barvec.^2,2));                        % bar lengths
barmid     = (GCOORD(bars(:,1),:) + GCOORD(bars(:,2),:))/2; % midpoint of each bar
switch SETTINGS.refinement
    case 'regular'
        L0 = SETTINGS.h0*ones(size(barmid,1),1);
    case 'guide_mesh'
        L0 = bar_L0_guide(barmid,GUIDE_MESH);
end

end % END OF FUNCTION bar_length