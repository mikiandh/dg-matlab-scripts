function groupAnimatedLines(axes)
% Assumming that all axes given contain (only) animated lines, plot as many
% figures as different lines there are in the first axis, and group in each
% figure all animated lines with a common local index (i.e. first with
% first, and so on).
%
% Setup:
N = numel(axes(1).Children);
M = numel(axes);
colors = distinguishable_colors(M,{'k','w'});

% Loop over lines per axis:
for i = 1:N
    figure()
    hold all
    % Loop over axes:
    for j = 1:M
        line = axes(j).Children(i);
        [x,y] = getpoints(line);
        plot(x,y,'DisplayName',sprintf('%s; figure %d, line %d',line.DisplayName,j,i),...
            'Color',colors(j,:),'LineStyle',line.LineStyle,'Marker',line.Marker)
    end
    legend('Location','Best')
end
end
        