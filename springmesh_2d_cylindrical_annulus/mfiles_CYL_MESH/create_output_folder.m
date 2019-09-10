function SETTINGS = create_output_folder(SETTINGS)
% Usage: SETTINGS = create_output_folder(SETTINGS)
%
% Purpose: Checks if output folder "SETTINGS.outdir"exists; 
%          creates it if necessary
%
% Input:
%   SETTINGS : [structure] : model parameters
%
% Output:
%   SETTINGS : [structure] : model parameters
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

try
    a=ls(SETTINGS.outdir);
catch %#ok<CTCH>
    a = [];
end
if isempty(a)
    try
        mkdir(SETTINGS.outdir);
    catch %#ok<CTCH>
        error(' Could not create output directory. Check if parent directory exists.');
    end
    fprintf(1,' Created output folder\n "%s" \n',SETTINGS.outdir);
else
    fprintf(1,' Output folder\n "%s" \n exists and data will be overwritten!!! \n',...
            SETTINGS.outdir);
end

if SETTINGS.outdir(end)~='/'
    SETTINGS.outdir = [SETTINGS.outdir '/'];
end

if isunix
    cmd = ['cp -p "./triangle_unix" "' SETTINGS.outdir '/triangle"'];
    system(cmd);
end

end % END OF FUNCTION create_output_folder