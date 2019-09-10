function display_add_reject_nodes(ADDREJECT)
% Usage: display_add_reject_nodes(ADDREJECT)
% 
% Purpose: Displays added and rejected nodes in structure "ADDREJECT"
%           in Matlab's terminal.
%
% Input:
%   ADDREJECT : [structure] : number of added and rejected nodes
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
% JMT May 2015: Addapted to MESH_2D_SPRING_CYL
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

varnames = fieldnames(ADDREJECT);
nvar     = length(varnames);
fprintf('\n Nodes added and rejected...\n');
fprintf(' --------------------------------\n');
fprintf(' Region  |  added  |  rejected  |\n');
for i=1:nvar
    varname = varnames{i};
    l       = 8 - length(varname);
    m       = 7 - length(num2str(ADDREJECT.(varname)(1)));
    n       = 10 - length(num2str(ADDREJECT.(varname)(2)));
    fprintf(' %s%s| ',varname,repmat(' ',1,l));
    fprintf(' %i%s|  %i%s|\n',ADDREJECT.(varname)(1),repmat(' ',1,m),ADDREJECT.(varname)(2),repmat(' ',1,n));
end
fprintf(' --------------------------------');

end % END OF FUNCTION display_add_reject_nodes
