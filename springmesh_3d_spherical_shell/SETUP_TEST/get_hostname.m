function hostname = get_hostname()
% Usage: hostname = get_hostname()
% 
% Purpose: Returns name of computer (i.e. host name)
%
% Input:
%   none
%
% Output:
%   hostname : [char] : name of computer
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2013
%

if strfind(pwd,'/Users/joha/')==1
    hostname = 'snaefell';
elseif strfind(pwd,'/Users/jorge/Dropbox/GEOMAR/springmesh_3d')==1
    hostname = 'jorgeMAC';
elseif strfind(pwd,'C:\Users\Tara\Documents\Dropbox\GEOMAR\springmesh_3d')==1
    hostname = 'jorgeLAPTOP';
elseif strfind(pwd,'/Users/jpm/Desktop/m2tri_cyl(trunk)')==1
    hostname = 'jasonLAPTOP';
elseif strfind(pwd,'/home/mat1/Jorge/springmesh_3d')==1
    hostname = 'clusterRHUL';
else
    tmp=pctconfig;
    hostname=tmp.hostname;
end

end