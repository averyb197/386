function astyle(colormap_name)
    if nargin < 1
        colormap_name = 'jet';  
    end
    
    numColors = 7;  % Number of colors to extract for line/marker cycling
    vibrantColors = feval(colormap_name, numColors);  

    set(0, 'DefaultFigureColor', [0, 0, 0]);          
    set(0, 'DefaultAxesColor', [0, 0, 0]);            
    set(0, 'DefaultAxesXColor', [1, 1, 1]);          
    set(0, 'DefaultAxesYColor', [1, 1, 1]);          
    set(0, 'DefaultAxesZColor', [1, 1, 1]);          
    set(0, 'DefaultAxesGridColor', [1, 1, 1]);       
    set(0, 'DefaultAxesGridAlpha', 0.2);             
    set(0, 'DefaultTextColor', [1, 1, 1]);
    set(0, 'DefaultAxesFontSize', 12);                
    set(0, 'DefaultAxesFontWeight', 'bold');         

    set(0, 'DefaultAxesXGrid', 'on');
    set(0, 'DefaultAxesYGrid', 'on');
    set(0, 'DefaultAxesZGrid', 'on');

    set(0, 'DefaultAxesColorOrder', vibrantColors);
    set(0, 'DefaultLineLineWidth', 1.5);              
    set(0, 'DefaultLineMarkerSize', 8);              
    colormap(colormap_name);
    message = sprintf("Hi Norberto, now your graphs will look cool. You're welcome. \n" + ...
        "Also, thank you for a great semester,\n this has by far been my favorite class of all time. ");
    disp(message)
end
