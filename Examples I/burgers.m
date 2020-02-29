clc
clear
%close all
%path(pathdef)

% This script solves the (inviscid) Burgers' equation.

%% Dependencies
addpath('../Extra')
addpath('../Discretization')
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
Ne = 20; % number of elements (!-> for 1 element, use FIXED time-step size)
p = 2; % degree of the approximation space (per element)
L = [-3 3]; % domain edges
tEnd = .5; % final simulation time
dt = [];
CFL = .2; % Courant number
iterSkip = 1;

%% Discretization
method = DG;

%% Physics
fun = @leveque1; % initial condition
eqn = Burgers(fun(L)); % PDE

%% Limiter
limiter = TVB('M',0);

%% Grid
mesh = Mesh(linspace(L(1),L(2),Ne+1),p,method,eqn);

%% Initial condition projection
%method.interpolate(mesh,limiter,fun);
method.project(mesh,limiter,fun);
%method.projectLumped(mesh,limiter,fun);

%% Time-integration
solver = SSP_RK3(0,tEnd,CFL,dt,limiter);
% norms0 = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation];
norms0 = [mesh.getSolutionMass mesh.getTVM];
tic
solver.launch(mesh,iterSkip,@(t,x) fun(x));
timeCPU = toc;
disp(['   ...done. (' num2str(timeCPU) ' s)'])

%% Postprocessing
% norms = [mesh.getSolutionMass,mesh.getSolutionNorm(1),mesh.getSolutionNorm,mesh.getTotalVariation];
% rows = {'Solution (mass)' 'Solution (L1)' 'Solution (L2)','Solution (TV)'};
norms = [mesh.getSolutionMass mesh.getTVM];
rows = {'Solution (mass)', 'Solution (TVM)'};
cols = {'Norm','Start','End','Increase'};
eqn.displayData(rows,cols,norms0,norms,norms-norms0)

%%% Plot the exact solution to Leveque1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(-3,3,1000);
y = @(x) -1 + (1+x/tEnd).*heaviside(x+tEnd) + (2-x/tEnd).*heaviside(x-2*tEnd);
hold all
h = get(gca,'Children');  
h(end).LineStyle = ':';
plot(x,y(x),'-.k')
fprintf('L2 norm of the error: %g \n',mesh.getErrorNorm(y,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function y = leveque1(x) % Leveque (expansion wave), 2002 (p. 231); L = [-3 3], tEnd = .5.
% NOTE: cancellation of fluxes may occur if -qL = qR
% NOTE: IGA breaks whenever N is even if shock is at x == 0
y = -1 + 3*heaviside(x);
end

function y = leveque2(x) % Adapted from Leveque, 2002 (p. 238); L = [-1 15], tEnd = 14.
y = 2 - 2*heaviside(x);
end

function y = leveque3(x) % Inspired by Leveque, 2002 (p. 223); L = [-8 8], tEnd = 6.
y = Functions.gauss(x,-4,4).*Functions.tones(x,[1 4 8],[1 2 3],-8,8);
end

function y = hesthaven(x) % Hesthaven & Warburton, 2008 (p. 141); L = [-1 1], tEnd = 0.8.
y = 2 - heaviside(x+.5);
end