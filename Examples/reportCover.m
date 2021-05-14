clc
clear

% This script generates a figure in A4 size and exports it to pdf; it could
% be used to generate a full-page cover image for a thesis.
%
% Requires 'export_fig' (https://github.com/altmany/export_fig) and 
% 'cbrewer' (https://uk.mathworks.com/matlabcentral/fileexchange/
% 34087-cbrewer-colorbrewer-schemes-for-matlab).

%% Page setup
X = 21.0;                  % A4 paper width
Y = 29.7;                  % A4 paper height
xMargin = .4;              % left/right margins from page borders
yMargin = .5;              % bottom/top margins from page borders
xSize = X - 2*xMargin;     % figure size on paper (widht & height)
ySize = Y - 2*yMargin;     % figure size on paper (widht & height)

%% Figure
hFig = figure;
basis = DGIGA(2,3,1);
cmap = brewermap(basis.basisCount,'Dark2');

% Optimal B-spline basis of 6 DoFs
ax1 = axes('Parent',hFig,'Units','normalized','Position',[.1 .2 .3 .2],...
    'ColorOrder',cmap,'NextPlot','replacechildren');
x = linspace(-1,1,1000);
plot(ax1,x,basis.sampleAt(x),'LineWidth',1);
ylabel('$N$','Interpreter','latex')
xlabel('$\xi$','Interpreter','latex')

% All its wavenumbers (Bloch waves)
ax2 = axes('Parent',hFig,'Units','normalized','Position',[.5 .1 .4 .4],...
    'ColorOrder',cmap,'NextPlot','replacechildren');
basis.displayModifiedWavenumbers;
set(ax2,...
    'XTick',[-pi 0 pi],'XTickLabel',{'-\pi','0','\pi'},...
    'YTick',[-pi 0 pi],'YTickLabel',{'-\pi','0','\pi'})
set(ax2.Children,'LineWidth',1)
% view(70,40)

%% Export
% Size on screen
set(hFig,'Units','centimeters','Position',[0 0 xSize ySize])

% Size on paper
set(hFig,'PaperUnits','centimeters')
set(hFig,'PaperSize',[X Y])
set(hFig,'PaperPosition',[xMargin yMargin xSize ySize])
set(hFig,'PaperOrientation','portrait')

% Background
set(hFig,'color','none');
set([ax1 ax2],{'color','SortMethod'},{'none','ChildOrder'});

export_fig('optimal.pdf','-nocrop','-transparent')