function MESH = plot_first_iter_reaching_different_worst_q_values(MESH,SETTINGS)
% Usage: MESH = plot_first_iter_reaching_different_worst_q_values(MESH,SETTINGS)
%
% Purpose:
%   Plot the 1st iteration where the worst quality factor reaches 0.20,
%   0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55 and 0.60 
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%   SETTINGS : [structure] : structure containing mesh settings
%
% Output:
%   MESH     : [structure] : structure containing the mesh
%
% JMT Jul 2016
%
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%--------------------------------------------------------------------------

if ~SETTINGS.show_figs
    set(gcf,'Visible','off');
end

% save 1st iter where worst_q reaches different values
MESH = save_first_iter_reaching_different_worst_q_values(MESH);

x         = MESH.iter_num;
x(x == 0) = NaN;
y         = MESH.worst_quality;
y(y == 0) = NaN;
plot(x,y,'r+')
axis([0 SETTINGS.itmax 0.2 0.65])
xlabel('number of iteration')
ylabel('quality factor')
title('First iteration in which worst q reaches different ranges of quality factor')

if SETTINGS.save_figs && (sum(~isnan(y)) == MESH.counter)
    prefix  = [SETTINGS.outdir '/'];
    suffix  = ['_' num2str_d(SETTINGS.iplot,5)];
    switch SETTINGS.fig_type
        case 'png'
            print('-dpng','-r150',[prefix 'Worst_q_evolution' suffix]);
        case 'fig'
            saveas(gcf,[prefix 'Worst_q_evolution' suffix],'fig');
        otherwise
            error(' Invalid value of "SETTINGS.fig_type".');
    end
    MESH.counter = MESH.counter + 1;
end

end % END OF FUNCTION plot_first_iter_reaching_different_worst_q_values

function MESH = save_first_iter_reaching_different_worst_q_values(MESH)
% Usage: MESH = save_first_iter_reaching_different_worst_q_values(MESH)
%
% Purpose:
%   Save the 1st iteration where the worst quality factor reaches 0.20,
%   0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55 and 0.60 
%
% Input:
%   MESH     : [structure] : structure containing the mesh
%
% Output:
%   MESH     : [structure] : structure containing the mesh
%
% JMT Jul 2016

worst_quality = MESH.worst_quality; % worst q once one of the values (0.20
                                    % 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                                    % 0.55 and 0.60) is reached
iter_num      = MESH.iter_num; % iteration in which worst_quality is computed
worst_q       = MESH.worst_q;
iter          = MESH.iter;
opt           = MESH.opt;

if sum(worst_quality)==0
    opt = '1st';
end

switch opt
    case '1st'
        if     worst_q >= 0.20 && worst_q < 0.25
            
            worst_quality(1)   = worst_q;
            iter_num(1)        = iter;
            opt                = '2nd';
            
        elseif worst_q >= 0.25 && worst_q < 0.30
            
            worst_quality(1)   = NaN;
            iter_num(1)        = NaN;
            worst_quality(2)   = worst_q;
            iter_num(2)        = iter;
            opt                = '3rd';
            
        elseif worst_q >= 0.30 && worst_q < 0.35
            
            worst_quality(1:2) = NaN;
            iter_num(1:2)      = NaN;
            worst_quality(3)   = worst_q;
            iter_num(3)        = iter;
            opt                = '4th';
            
        elseif worst_q >= 0.35 && worst_q < 0.40
            
            worst_quality(1:3) = NaN;
            iter_num(1:3)      = NaN;
            worst_quality(4)   = worst_q;
            iter_num(4)        = iter;
            opt                = '5th';
            
        elseif worst_q >= 0.40 && worst_q < 0.45
            
            worst_quality(1:4) = NaN;
            iter_num(1:4)      = NaN;
            worst_quality(5)   = worst_q;
            iter_num(5)        = iter;
            opt                = '6th';
            
        elseif worst_q >= 0.45 && worst_q < 0.50
            
            worst_quality(1:5) = NaN;
            iter_num(1:5)      = NaN;
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(1:6) = NaN;
            iter_num(1:6)      = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(1:7) = NaN;
            iter_num(1:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(1:8) = NaN;
            iter_num(1:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
            
        end
        
    case '2nd'
        if     worst_q >= 0.25 && worst_q < 0.30
            
            worst_quality(2)   = worst_q;
            iter_num(2)        = iter;
            opt                = '3rd';
            
        elseif worst_q >= 0.30 && worst_q < 0.35
            
            worst_quality(2)   = NaN;
            iter_num(2)        = NaN;
            worst_quality(3)   = worst_q;
            iter_num(3)        = iter;
            opt                = '4th';
            
        elseif worst_q >= 0.35 && worst_q < 0.40
            
            worst_quality(2:3) = NaN;
            iter_num(2:3)      = NaN;
            worst_quality(4)   = worst_q;
            iter_num(4)        = iter;
            opt                = '5th';
            
        elseif worst_q >= 0.40 && worst_q < 0.45    
            
            worst_quality(2:4) = NaN;
            iter_num(2:4)      = NaN;
            worst_quality(5)   = worst_q;
            iter_num(5)        = iter;
            opt                = '6th';
            
        elseif worst_q >= 0.45 && worst_q < 0.50
            
            worst_quality(2:5) = NaN;
            iter_num(2:5)      = NaN;
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(2:6) = NaN;
            iter_num(2:6)      = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(2:7) = NaN;
            iter_num(2:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(2:8) = NaN;
            iter_num(2:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
                
        end
        
    case '3rd'
        if     worst_q >= 0.30 && worst_q < 0.35
            
            worst_quality(3)   = worst_q;
            iter_num(3)        = iter;
            opt                = '4th';
            
        elseif worst_q >= 0.35 && worst_q < 0.40
            
            worst_quality(3)   = NaN;
            iter_num(3)        = NaN;
            worst_quality(4)   = worst_q;
            iter_num(4)        = iter;
            opt                = '5th';
            
        elseif worst_q >= 0.40 && worst_q < 0.45
            
            worst_quality(3:4) = NaN;
            iter_num(3:4)      = NaN;
            worst_quality(5)   = worst_q;
            iter_num(5)        = iter;
            opt                = '6th';
            
        elseif worst_q >= 0.45 && worst_q < 0.50
            
            worst_quality(3:5) = NaN;
            iter_num(3:5)      = NaN;
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(3:6) = NaN;
            iter_num(3:6)      = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(3:7) = NaN;
            iter_num(3:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(3:8) = NaN;
            iter_num(3:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
                
        end
        
    case '4th'
        if     worst_q >= 0.35 && worst_q < 0.40
            
            worst_quality(4)   = worst_q;
            iter_num(4)        = iter;
            opt                = '5th';
            
        elseif worst_q >= 0.40 && worst_q < 0.45
            
            worst_quality(4)   = NaN;
            iter_num(4)        = NaN;
            worst_quality(5)   = worst_q;
            iter_num(5)        = iter;
            opt                = '6th';
            
        elseif worst_q >= 0.45 && worst_q < 0.50
            
            worst_quality(4:5) = NaN;
            iter_num(4:5)      = NaN;
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(4:6) = NaN;
            iter_num(4:6)      = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(4:7) = NaN;
            iter_num(4:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(4:8) = NaN;
            iter_num(4:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
                
        end
        
    case '5th'
        if     worst_q >= 0.40 && worst_q < 0.45
            
            worst_quality(5)   = worst_q;
            iter_num(5)        = iter;
            opt                = '6th';
            
        elseif worst_q >= 0.45 && worst_q < 0.50
        
            worst_quality(5)   = NaN;
            iter_num(5)        = NaN;
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(5:6) = NaN;
            iter_num(5:6)      = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(5:7) = NaN;
            iter_num(5:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(5:8) = NaN;
            iter_num(5:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
                
        end
        
    case '6th'
        if     worst_q >= 0.45 && worst_q < 0.50
            
            worst_quality(6)   = worst_q;
            iter_num(6)        = iter;
            opt                = '7th';
            
        elseif worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(6)   = NaN;
            iter_num(6)        = NaN;
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(6:7) = NaN;
            iter_num(6:7)      = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
         elseif worst_q >= 0.60
            
            worst_quality(6:8) = NaN;
            iter_num(6:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
                
        end
        
    case '7th'
        if     worst_q >= 0.50 && worst_q < 0.55
            
            worst_quality(7)   = worst_q;
            iter_num(7)        = iter;
            opt                = '8th';
            
        elseif worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(7)   = NaN;
            iter_num(7)        = NaN;
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(7:8) = NaN;
            iter_num(7:8)      = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
            
        end
        
    case '8th'
        if     worst_q >= 0.55 && worst_q < 0.60
            
            worst_quality(8)   = worst_q;
            iter_num(8)        = iter;
            opt                = '9th';
            
        elseif worst_q >= 0.60
            
            worst_quality(8)   = NaN;
            iter_num(8)        = NaN;
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
            
        end
        
    case '9th'
        if     worst_q >= 0.60
            
            worst_quality(9)   = worst_q;
            iter_num(9)        = iter;
        end
end

MESH.worst_quality = worst_quality;
MESH.iter_num      = iter_num;
MESH.opt           = opt;

end % END OF SUBFUNCTION save_first_iter_reaching_different_worst_q_values