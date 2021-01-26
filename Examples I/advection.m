clc
clear
%close all

% This script solves the advection equation in 1D.

%% Parameters
L = [-1 1]; % domain edges

%% Initial condition collection
IC_linear = @(x) x;
IC_quadratic = @(x) 2*(x + x.^2);
IC_heaviside = @(x) heaviside(2/diff(L)*(x-L(1)) - 1);
IC_gauss = @(x) exp(-18*(.5*sqrt(2*pi))^2*(x-.5*(L(1)+L(2))).^2/(L(2)-L(1))^2);
IC_gaussgauss = @(x) IC_gauss(2*(x-L(1))+L(1)) + IC_gauss(2*(x-.5*sum(L))+L(1));
IC_sine = @(x) 1+.5*sin(2*pi*(x)/diff(L));
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

%% Discretization
mesh = Mesh(DGSEM(5),L,Periodic(2),100);

%% Solver
solver = SSP_RK4_10(Advection,[0 2],...
    'limiters',[WENO Limiter Limiter],...
    'exactSolution',@(t,x) IC_jiangShu(x),...
    'iterSkip',25);
solver.courantNumber = .1*solver.optimizeCFL(mesh.bases);

%% Initial condition
solver.initialize(mesh)

%% Time-integration
solver.launch(mesh)