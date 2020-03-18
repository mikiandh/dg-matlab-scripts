clc
clear
%close all
%path(pathdef)

% This script solves the advection equation in 1D.

%% Dependencies
addpath('../Extra')
addpath('../Discretization')
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
L = [-1 1]; % domain edges

%% Initial condition collection
IC_linear = @(x) x;
IC_quadratic = @(x) 2*(x + x.^2);
IC_heaviside = @(x) heaviside(2/diff(L)*(x-L(1)) - 1);
IC_gauss = @(x) exp(-18*(.5*sqrt(2*pi))^2*(x-.5*(L(1)+L(2))).^2/(L(2)-L(1))^2);
IC_gaussgauss = @(x) IC_gauss(2*(x-L(1))+L(1)) + IC_gauss(2*(x-.5*sum(L))+L(1));
IC_sine = @(x) sin(2*pi*(x)/diff(L));
IC_jump = @(x) heaviside(x-2/3*L(1)-1/3*L(2)) - heaviside(x-1/3*L(1)-2/3*L(2));
IC_opposedJumps = @(x) IC_jump(2*(x-L(1))+L(1)) - IC_jump(2*(x-.5*sum(L))+L(1));
IC_packet = @(x) wavePacket(x,diff(L)/2,[2],[.5]'.*IC_gauss(x),0.5); %#ok
IC_combined = @(x) IC_gauss(2*(x-L(1))+L(1)) + IC_jump(2*(x-.5*sum(L))+L(1));
IC_jiangShu = @(x) JiangShu(2/diff(L)*(x-L(1)) - 1);
IC_randomSignal = @(x) randomSignal(x,0,3,L,0);
IC_jaeschkeHump = @(x) .25*(1+cos(pi*min(1,2*sqrt((x-.5).^2))));
IC_jaeschkeSquare = @(x) heaviside(x-0.25) - heaviside(x-0.75);
IC_leveque = @(x) 2 - 2*heaviside(x);
IC_p3d0 = @(x) (x+x.^2+x.^3).*(1 - heaviside(x)) + (1+x+x.^2+x.^3).*(heaviside(x));
IC_p3d1 = @(x) (x+x.^2+x.^3).*(1 - heaviside(x)) + (2*x+x.^2+x.^3).*(heaviside(x));
IC_p3d2 = @(x) (x+x.^2+x.^3).*(1 - heaviside(x)) + (x+2*x.^2+x.^3).*(heaviside(x));
IC_p3d3 = @(x) (x+x.^2+x.^3).*(1 - heaviside(x)) + (x+x.^2+2*x.^3).*(heaviside(x));
IC_p3d4 = @(x) (x+x.^2+x.^3);

%% Discretization (equation + solution + domain)
xEdge = linspace(L(1),L(2),10+1); % element end-points
mesh = Mesh(xEdge,2,DGIGA_AFC(10));

%% Solver
solver = SSP_RK3(0,0,Advection(1,[]),...
    'courantNumber',.01,...
    'limiter',AFC_2010('Sensor',APTVD),...
    'exactSolution',@(t,x) IC_combined(x),'replotIters',1);

%% Initial condition
solver.initialize(mesh)
% norms0 = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation,mesh.getErrorNorm(FUN)];
norms0 = mesh.getSolutionMass;

%% Time-integration
solver.launch(mesh);

%% Postprocessing
% norms = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation,mesh.getErrorNorm(FUN)];
% rows = {'Solution (mass)' 'Solution (L1)' 'Solution (L2)','Solution (TV)','Error (L2)'};
norms = mesh.getSolutionMass;
rows = {'Solution (mass)'};
cols = {'Norm','Start','End','Increase'};
solver.physics.displayData(rows,cols,norms0,norms,norms-norms0)