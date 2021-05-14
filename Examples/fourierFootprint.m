clc
clear

% This script is meant to allow a visual comparison between the Fourier 
% footprint of various spatial discretizations in relation to the stability 
% regions of a single time discretization.

%% Spatial discretizations
discs = {
    DGSEM(5)
    FR({'eta',0.058},5)
    DGIGA(2,3,1)
    };

%% Temporal discretization
time = SSP_RK3;

%% Extra parameters
% Upwinding ratio for Riemann fluxes, from -1 (downwind) to 1 (upwind):
upwind = 1;
% Number of patches (only affects resolution; should be even):
Nx = 32;

%% Fourier footprints
I = size(discs,1);
z = cell(I,1);
CFL = ones(I,1);
for i = 1:I
    % Get a Fourier footprint:
    z{i} = complex(pi*discs{i}.basisCount*linspace(0,1,Nx));
    z{i} = reshape(discs{i}.getFourierFootprint('upwind',upwind,'wavenumbers',z{i}),1,[]);
    % Scale Fourier footprint with its critical Courant number:
    CFL(i) = time.optimizeCFL(discs{i},'upwind',upwind);
    z{i} = z{i}*CFL(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % z{i} = z{i}*0.0525;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recover negative wavenumbers (aesthetic reasons only):
    z{i} = horzcat(flip((z{i}').',2),z{i});
end

%% Visualization
figure('Renderer', 'painters', 'Position', [100 100 800 800])
c = distinguishable_colors(I,{'w','k',[.8 1 1]});
% Stability region:
[x,y] = meshgrid(linspace(-16,2,250),linspace(-8,8,250));
contourf(x,y,abs(time.amplificationFactorFun(x+1i*y)),[0 1]);
colormap([.8 1 1; 1 1 1]);
hold all
% Coordinate axes:
plot([x(1) x(end)],[0 0],'k--')
plot([0 0],[y(1) y(end)],'k--')
% Fourier footprints:
h = zeros(1,I);
for i = 1:I
    h(i) = plot(real(z{i}),imag(z{i}),...
        '.','Color',c(i,:),'DisplayName',discs{i}.getName);
end
hold off
setFancyPlotSettings3
axis equal
xlabel('$$\Re(\tilde{z})$$','Interpreter','latex')
ylabel('$$\Im(\tilde{z})$$','Interpreter','latex')
if isempty(time)
    title(strrep(char(G),'@(z)',''))
else
    title(class(time),'Interpreter','none')
end
legend(h,'Location','SouthWest')