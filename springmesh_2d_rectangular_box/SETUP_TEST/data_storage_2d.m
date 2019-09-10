function path2data = data_storage_2d()
% Usage: path2data = data_storage_2d()
% 
% Purpose: Sets path to data storage location on current computer.
%          (Edit this file to add new computers)
%
% Input:
%   none
%
% Output:
%   path2data : [char] : complete path to data storage
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2012

hostname = get_hostname();
switch hostname
    case 'b3pc17'
%         path2data = 'E:/Junk Yard/';
        path2data = 'D:/Projects/';
    case 'B3PC12'
        path2data = './OUTPUT/';
    case 'b3pc14'
        path2data = 'E:/Projects/';
    case 'b3pc16'
%         path2data = 'J:/jhasenclever/';
        path2data = 'Y:/';
    case 'b3pc9'
        path2data = 'E:/jhasenclever/';
%         path2data = 'Z:/';
    case 'Maupiti2'
        path2data = 'F:/Projects/';
    case 'snaefell'
        path2data = '/Users/joha/TMP/';
    case 'fuego'
        path2data = '/Data_Fu/jhasenclever/';
    case 'colima'
        path2data = '/Data_Co/jhasenclever/';
    case 'jorgeMAC'
        path2data = '/Users/jorge/Tests/';
    case 'jorgeLAPTOP'
        path2data = 'C:\Users\Tara\Documents\PhD\Tests\';
    case 'jasonLAPTOP'
        path2data = '/Users/jpm/Desktop/m2tri_cyl(trunk)';
    otherwise
        if strcmp(hostname(1:5),'node0') || strcmp(hostname,'master')
            path2data = './OUTPUT/';
        else
            error(' No location for output data defined for this computer.\n You have to edit file "data_storage_2d.m"!');
        end
end

% '\' can cause problems on unix system; however matlab can handle '/' on
% all architectures
path2data(path2data=='\') = '/';
% make sure there's a slash at the end of the path
if path2data(end)~='/'
    path2data = [path2data '/'];
end

end % END OF FUNCTION data_storage_2d