% this is a ruotine for output Matlab figures.
% cut the marginal white space.
% created on Sep. 27,2015; 

% this routine should be put between "figure" and "plot"

% width = 5.8; height = 6; is the "square"
% width = 7;   height = 6; is the rantangular, 4:3


function My_Figure( width, height )
    
    load MyColors;
    set(gca, 'LooseInset', get(gca, 'TightInset')); 
    %set(gca,'fontsize',24); % my default is 28
    
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [width-0.5 height-0.5]);
    set(gcf, 'PaperPosition', [0 0 width height]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'renderer', 'painters');
    
    %xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);
    %setappdata(0, 'DefaultAxesXLabelFontSize', 36); % no use
end

