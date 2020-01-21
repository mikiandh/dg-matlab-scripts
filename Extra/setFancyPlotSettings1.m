%% Plot makeup settings
% Tweaks the plot settings so that thwe result looks better. Version 1.
function setFancyPlotSettings1()
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'XGrid'       , 'on'      , ...
    'XMinorGrid'  , 'on'     , ...
    'YGrid'       , 'on'      , ...
    'YMinorGrid'  , 'off'     , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1);
set(gcf, 'Color', 'w');
end