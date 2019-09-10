function [pfix,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS)
% Usage: [pfix,pbnd1,pbnd2] = boundary_nodes_no_guide(SETTINGS)
%
% Purpose:
%   Generation of boundary nodes for a 2D cylindrical annulus regular mesh
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   pfix     : [matrix]    : coord of fixed nodes
%   pbnd1    : [matrix]    : coord of boundary nodes for inner boundary
%   pbnd2    : [matrix]    : coord of boundary nodes for outer boundary
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

%==========================================================================
% LOAD VARIABLES
%==========================================================================
r_int = SETTINGS.r_int;
r_ext = SETTINGS.r_ext;
h0    = SETTINGS.h0;

%==========================================================================
% BOUNDARY 1 (INNER BOUNDARY)
%==========================================================================
if r_int > 0
    pbnd1         = regular_bnd_nodes(r_int,h0);
    [pfix1,pbnd1] = create_pfix(pbnd1,r_int);
elseif r_int == 0
    pfix1 = [];
    pbnd1 = [];
else
    error('r_int must be >= 0')
end

%==========================================================================
% BOUNDARY 2 (OUTER BOUNDARY)
%==========================================================================
pbnd2         = regular_bnd_nodes(r_ext,h0);
[pfix2,pbnd2] = create_pfix(pbnd2,r_ext);
% pfix2(2,:) = [ r_ext 0];
% pfix2(4,:) = [-r_ext 0];
%==========================================================================
% OUPUT DATA
%==========================================================================
pfix           = [pfix1; pfix2];

if strcmp(SETTINGS.mesh,'axisym')
    if SETTINGS.r_int == 0;
        y          = linspace(r_ext,-r_ext,round(2*r_ext/h0))';
        y([1,end]) = [];
        x          = zeros(size(y,1),1);
        pfix(3,:)  = [0 -r_ext];
        pfix       = [pfix; [x y]];
        pfix(pfix(:,1) < 0,:)    = []; % remove pfix points with x < 0
        pbnd2(pbnd2(:,1) <= 0,:) = []; % remove pbnd2 points with x <= 0
    else
        y1          = linspace(r_ext,r_int,round((r_ext-r_int)/h0))';
        y1([1,end]) = [];
        x1          = zeros(size(y1,1),1);
        y2          = linspace(-r_int,-r_ext,round((r_ext-r_int)/h0))';
        y2([1,end]) = [];
        x2          = zeros(size(y2,1),1);
        pfix(3,:)   = [0 -r_int];
        pfix(7,:)   = [0 -r_ext];
        pfix        = [pfix; [x1 y1]; [x2 y2]];
        pfix(pfix(:,1) < 0,:)    = []; % remove pfix points with x < 0
        pbnd1(pbnd1(:,1) <= 0,:) = []; % remove pbnd1 points with x <= 0
        pbnd2(pbnd2(:,1) <= 0,:) = []; % remove pbnd2 points with x <= 0
    end
end

end % END OF FUNCTION boundary_nodes_no_guide