clc
clear
%close all
%path(pathdef)

% This script solves the (inviscid) Burgers' equation.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
Ne = 1; % number of elements (!-> for 1 element, use FIXED time-step size)
p = 1; % degree of the approximation space (per element)
L = [-1 9]; % domain edges
tEnd = 4; % final simulation time
dt = .25;
CFL = .1; % Courant number
iterSkip = 1;

%% Discretization
method = DGIGA_AFC(20);

%% Grid
mesh = Mesh(linspace(L(1),L(2),Ne+1),p,method);

%% Physics
FUN = @leveque2; % initial condition
eqn = Burgers(FUN(L)); % PDE

%% Limiter
limiter = [];
%limiter = Limiter.TVBM(eqn,0);
%limiter = Limiter.TVBM(eqn,5);
%limiter = Limiter.Biswas(eqn);
%limiter = Limiter.Burbeau(eqn);
%limiter = Limiter.Krivodonova(eqn);
%limiter = Limiter.Wang(eqn); % best in this case (close call)
if isa(method,'DGIGA_AFC')
    %%%limiter = Limiter.AFC(eqn);
end

%% Initial condition projection
% method.interpolate(mesh,limiter,FUN);
method.project(mesh,limiter,FUN);

%% Time-integration
timeIntegrator = SSP_RK1(0,tEnd,eqn,limiter,CFL,dt);
% norms0 = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation];
norms0 = mesh.getSolutionMass;
tic
timeIntegrator.launch(mesh,iterSkip,@(t,x) FUN(x));
disp(['   ...done. (' num2str(toc) ' s)'])

%% Postprocessing
% norms = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation];
% rows = {'Solution (mass)' 'Solution (L1)' 'Solution (L2)','Solution (TV)'};
norms = mesh.getSolutionMass;
rows = {'Solution (mass)'};
cols = {'Norm','Start','End','Increase'};
eqn.displayData(rows,cols,norms0,norms,norms-norms0)

%% Initial conditions
function y = riemann1A(x) % right-going shock
y = 1 - 0.8*heaviside(x);
end

function y = riemann1B(x) % right-going expansion
y = 0.2 + 0.8*heaviside(x);
end

function y = riemann2A(x) % left-going expansion
y = -1 + 0.8*heaviside(x);
end

function y = riemann2B(x) % left-going shock
y = 1 - 1.2*heaviside(x);
end

function y = riemann3A(x) % shock-shock collision (right-going wins)
y = 1 - 1.2*heaviside(x);
end

function y = riemann3B(x) % shock-shock collision (left-going wins)
y = 0.2 - heaviside(x);
end

function y = riemann3C(x) % shock-shock collision (standing shock)
y = 1 - 2*heaviside(x);
end

function y = riemann4A(x) % transonic expansion (right-going)
y = -0.2 + heaviside(x);
end

function y = riemann4B(x) % transonic expansion (left-going)
y = -0.8 + heaviside(x);
end

function y = gaussIC(x)
y = Functions.gauss(x,-1,1);
end

function y = jumpIC(x)
y = Functions.jump(x,-1,1);
end

function y = halfJumpIC(x)
y = x.*(heaviside(x) - heaviside(x-2));
end

function y = sineIC(x)
y = 0.5*sin(9*x*pi) + 0.5;
end

function y = combinedIC(x) % perfect initial condition; L = [-1,2], tEnd = 2
y = Functions.gauss(x,-1,0) + Functions.jump(x,0,1);
end

function y = toroIC(x) % from Toro, 2009 (p. 196); L = [0 1.5], tEnd = .5
y = -.5 + 1.5*heaviside(x-.5) - 1*heaviside(x-1);
end

function y = leveque1(x) % Leveque (expansion wave), 2002 (p. 231); L = [-3 3], tEnd = 1.
y = -1 + 3*heaviside(x);
end

function y = leveque2(x) % Adapted from Leveque, 2002 (p. 238); L = [-1 15], tEnd = 14.
y = 2 - 2*heaviside(x);
end

function y = leveque3(x) % Inspired by Leveque, 2002 (p. 223); L = [-8 8], tEnd = 6.
y = Functions.gauss(x,-4,4).*Functions.tones(x,[1 4 8],[1 2 3],-8,8);
end