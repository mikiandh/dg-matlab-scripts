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
Ne = 20; % number of elements
p = 11; % degree of the approximation space (per element)
L = [-.5 .5]; % domain edges
tEnd = 1; % final simulation time
dt = []; % time-step size (overrides CFL)
CFL = .01; % Courant number
iterSkip = 10;

%% Initial condition collection
IC_linear = @(x) x;
IC_quadratic = @(x) 2*(x + x.^2);
IC_heaviside = @(x) heaviside(2/diff(L)*(x-L(1)) - 1);
IC_gauss = @(x) exp(-18*(.5*sqrt(2*pi))^2*(x-.5*(L(1)+L(2))).^2/(L(2)-L(1))^2);
IC_gaussgauss = @(x) IC_gauss(2*(x-L(1))+L(1)) + IC_gauss(2*(x-.5*sum(L))+L(1));
IC_sine = @(x) 1+sin(2*pi*(x)/diff(L));
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
IC_p3d3 = @(x) (x+x.^2+x.^3).*(1 - heaviside(x)) + (x+x.^2+2*x.^3)  .*(heaviside(x));
IC_p3d4 = @(x) (x+x.^2+x.^3);

%% Physics
FUN = IC_jump; % initial condition
eqn = Advection(1,[]); % PDE + BCs

%% Discretization
method = DGSEM;

%% Limiter
limiter = TVB('M',0);

%% Grid
xEdge = linspace(L(1),L(2),Ne+1); % element end-points
mesh = Mesh(xEdge,p,method,eqn);

%% Initial condition projection
%method.interpolate(mesh,limiter,FUN);
method.project(mesh,limiter,FUN);
%method.projectLumped(mesh,limiter,FUN);
%method.project_Matthias(mesh,limiter,FUN);

%% Time-integration
solver = SSP_RK3(0,tEnd,CFL,dt,limiter);
%norms0 = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation,mesh.getErrorNorm(FUN)];
norms0 = mesh.getSolutionMass;
tic
solver.launch(mesh,iterSkip,@(t,x) FUN(x));
fprintf(1,'...done. (%g s)\n',toc)

%% Postprocessing
% norms = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation,mesh.getErrorNorm(FUN)];
% rows = {'Solution (mass)' 'Solution (L1)' 'Solution (L2)','Solution (TV)','Error (L2)'};
norms = mesh.getSolutionMass;
rows = {'Solution (mass)'};
cols = {'Norm','Start','End','Increase'};
eqn.displayData(rows,cols,norms0,norms,norms-norms0)