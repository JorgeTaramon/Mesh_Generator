function [pfix,pbnd,segBC,pseg] = boundary_nodes_no_guide(SETTINGS)
% Usage: [pfix,pbnd,segBC,pseg] = boundary_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generation of boundary nodes for a 2D rectangular regular mesh
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   pfix     : [matrix]    : coord of fixed nodes
%   pbnd     : [matrix]    : coord of boundary nodes
%   segBC    : [matrix]    : info of boundary segments
%   pseg     : [vector]    : pointID for bnd nodes
%
% JMT Jun 2017
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
depth      = SETTINGS.depth;         % domain depth (km)
length     = SETTINGS.length;        % domain length (km)
x0         = SETTINGS.x0;            % point around which the domain is created
z0         = SETTINGS.z0;            % point around which the domain is created
h0         = SETTINGS.h0;            % desired spring length (km) for a regular mesh

%==========================================================================
% DEFINE PFIX
%==========================================================================
pfix = [x0 - length/2  z0 - depth ; ...
        x0 + length/2  z0 - depth ; ...
        x0 + length/2  z0 ; ...
        x0 - length/2  z0 ]; % vertex points COUNTER-CLOCKWISE
nfix = size(pfix,1);
nseg = nfix;

%==========================================================================
% DEFINE EACH BOUNDARY SEGMENT
%==========================================================================
seg_connect    = [(1:nseg)',[(2:nseg)';1]];  % This is the connectivity matrix for the fixed points
segBC          = zeros(nseg,6);              % pre-allocate memory
segBC(:,[1 2]) = pfix(seg_connect(:,1),:);   % p0 = (x0,y0) for seg_start
segBC(:,[3 4]) = pfix(seg_connect(:,2),:);   % p1 = (x1,y1) for seg_end
segBC(:,[5 6]) = GetLineSlope(segBC(:,[1 2]),segBC(:,[3 4]),nseg); % slope_info = [a,b]

%==========================================================================
% COMPUTE PBND
%==========================================================================
pbnd = [];
pseg = [];
for i = 1:nseg
    pbnd_this_segment = regular_bnd_nodes(segBC,h0,i);
    pbnd              = [pbnd; pbnd_this_segment];
    pseg              = [pseg i*ones(1,size(pbnd_this_segment,1))];
end

end % END OF FUNCTION boundary_nodes_no_guide