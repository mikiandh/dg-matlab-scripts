%% Plot makeup settings
% Tweaks the plot settings so that thwe result looks better. Version 2.
function setFancyPlotSettings2()
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'on'      , ...
    'XMinorGrid'  , 'on'     , ...
    'YGrid'       , 'on'      , ...
    'YMinorGrid'  , 'on'     , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1);
set(gcf, 'Color', 'w');
end