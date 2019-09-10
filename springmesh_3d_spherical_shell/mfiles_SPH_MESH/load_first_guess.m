function MESH = load_first_guess(SETTINGS)
% Usage: MESH = load_first_guess(SETTINGS)
%
% Purpose:
%   Load a first guess mesh
%
% Input:
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   MESH     : [structure] : structure containing 1st guess for the mesh
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

filename                = [pwd '/' SETTINGS.first_guess_file];
filename(filename=='\') = '/';
try
    tmp = load(filename);
catch
    error('Cannot find mesh file %s',filename);
end
if ~isfield(tmp,'MESH')
    error('Meshfile %s does not contain structure MESH.',filename);
end
MESH = tmp.MESH; clear tmp

if SETTINGS.save_figs || SETTINGS.show_figs
    plot_first_guess_mesh(MESH,SETTINGS)
end

end % END OF FUNCTION load_first_guess