function display_progress(iter,t0)
% Usage: display_progress(iter,t0)
% 
% Purpose: 
%   Display current model iteration and run time in terminal.
%
% Input:
%   iter : [scalar]    : current iteration
%   t0   : [rowvector] : Matlab's clock command from start of code
%
% Output:
%   none (Output only in Matlab terminal)
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
% JMT May 2016: working for MESH_2D_SPRING
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

fprintf('\n................  Iteration %4i  ................\n',iter);

fprintf('Code running since %2i d   %2i h   %2i min  %5.2f sec\n',...
        return_runtime(t0));

end % END OF FUNCTION display_progress

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function runtime = return_runtime(t0)

s = etime(clock,t0);

% For debugging
%     d = 1;
%     h = 4;
%     m = 32;
%     s = 17.3;
%     s = s + 60*m + 60*60*h + 60*60*24*d;
% The function should return [d h m s].

m = floor(s/60);
h = floor(m/60);
d = floor(h/24);
h = h-24*d;
m = m-60*h-60*24*d;
s = s-60*m-60*60*h-60*60*24*d;
runtime = [d h m s];
    
end % END OF SUBFUNCTION return_runtime