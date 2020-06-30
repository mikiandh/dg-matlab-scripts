clc
clear
%close all
%path(pathdef)

% This script shows the Fourier footprint of various spatial
% discretizations in relation to the stability region of a single time
% discretization.

%% Dependencies
addpath('../../Extra')
addpath('../../Basis')
addpath('../../Solver')

%% Spatial discretizations
discs = {...
    DGIGA(2,2,0)
    };

%% Temporal discretization
time = SSP_RK3;

%% Extra parameters
% Upwinding ratio for Riemann fluxes, from -1 (downwind) to 1 (upwind):
upwind = 1;
% Number of patches (only affects resolution):
Nx = 32;
% Temporal amplification factor function:
G = time.amplificationFactorFun;
% G = @(z) 1./(1 - z); % overwrite with implict RK1
% Courant number seed:
CFL = 1;
% Plot footprint as points instead of curves? (much more efficient)
flagPoints = true;

%% Fourier footprints
I = size(discs,1);
z = cell(I,1);
CFL = CFL*ones(I,1);
for i = 1:I
    % Modified wavenumber analysis:
    z{i} = reshape(...
        -1i*MWA_eigen_full(Nx,discs{i}.degree,discs{i},upwind),1,[]);
    % Scale Fourier footprint with its critical Courant number:
    CFL(i) = optimizeCFL(CFL(i),z{i},G);
    z{i} = z{i}*CFL(i);
    % Recover negative wavenumbers (aesthetic reasons only):
    z{i} = horzcat(flip((z{i}').',2),z{i});
    if ~flagPoints
        % Split:
        z{i} = vertcat(real(z{i}),imag(z{i}));
        % Reorder:
        [z{i}(1,:), z{i}(2,:)] = points2contour(z{i}(1,:),z{i}(2,:),1,'cw');
        % Stitch the ends together:
        z{i} = horzcat(z{i},z{i}(:,1));
    end
end

%% Visualization
figure('Renderer', 'painters', 'Position', [100 100 800 800])
c = distinguishable_colors(I,{'w','k',[.8 1 1]});
% Stability region:
[x,y] = meshgrid(linspace(-16,2,250),linspace(-8,8,250));
contourf(x,y,abs(G(x+1i*y)),[0 1]);
colormap([.8 1 1; 1 1 1]);
hold all
% Coordinate axes:
plot([x(1) x(end)],[0 0],'k--')
plot([0 0],[y(1) y(end)],'k--')
% Fourier footprints:
h = zeros(1,I);
for i = 1:I
    if isa(discs{i},'FR')
        name = sprintf('%s(%s), p = %d; \\varsigma = %.4f',...
            class(discs{i}),num2str(discs{i}.param),discs{i}.degree,CFL(i));
    elseif discs{i}.isHybrid
        name = sprintf('%s, p = %d, Ns = %d, C^{%d}; \\varsigma = %.4f',...
            class(discs{i}),discs{i}.degree,discs{i}.nonzeroSpanCount,...
            discs{i}.smoothness,CFL(i));
    else
        name = sprintf('%s, p = %d; \\varsigma = %.4f',...
            class(discs{i}),discs{i}.degree,CFL(i));
    end
    if flagPoints
        h(i) = plot(real(z{i}),imag(z{i}),...
            '.','Color',c(i,:),'DisplayName',name);
    else
        h(i) = plot(z{i}(1,:),z{i}(2,:),...
            '-','Color',c(i,:),'DisplayName',name); %#ok<UNRCH>
    end
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

%% Nonlinear constrained optimization subroutine
function [CFL,exitFlag] = optimizeCFL(CFL,theta,G)
% This function sets up a constrained optimization problem, and solves it 
% via fmincon. The goal is to determine the CFL coefficient that maximizes
% the L1 norm of the amplification factors - all eigenmodes and wavenumbers
% considered.
%
% Nonlinear function to minimize:
    function norm = fun(CFL)
        norm = - sum(abs(G(theta*CFL)));
    end
% Nonlinear constraint:
    function [c,ceq] = nonlcon(CFL)
        c = max(abs(G(theta*CFL))) - 1;
        ceq = [];
    end
% Solve the problem:
[CFL,~,exitFlag] = fmincon(@fun,CFL,[],[],[],[],0,[],@nonlcon,...
    optimoptions('fmincon','Display','off'));
end