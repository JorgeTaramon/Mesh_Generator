function GUIDE_MESH = load_guide_mesh(SETTINGS)
% Usage: GUIDE_MESH = load_guide_mesh(SETTINGS)
%
% Purpose:
%   Load a guide mesh
%
% Input:
%   SETTINGS   : [structure] : structure containing mesh settings
%
% Output:
%   GUIDE_MESH : [structure] : structure containing guide mesh
%
% JMT May 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

filename                = [pwd '/' SETTINGS.guide_mesh_file];
filename(filename=='\') = '/';
try
    tmp = load(filename);
catch
    error('Cannot find mesh file %s',filename);
end
if ~isfield(tmp,'GUIDE_MESH')
    error('Meshfile %s does not contain structure GUIDE_MESH.',filename);
end
GUIDE_MESH = tmp.GUIDE_MESH; clear tmp

if SETTINGS.save_figs || SETTINGS.show_figs
    plot_guide_mesh(GUIDE_MESH,SETTINGS)
end

end % END OF FUNCTION load_guide_mesh
    