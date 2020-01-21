clc
clear
%close all
%path(pathdef)

% This script solves the wave equation in 1D. Uses the DGSEM classes.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
Ne = 15; % number of elements (!-> for 1 element, use FIXED time-step size)
p = 4; % degree of the approximation space (per element)
L = [0 1]; % domain edges
tEnd = 4.0; % final simulation time
dt = [];
CFL = .05; % Courant number
iterSkip = 120;
    
%% Physics
eqn = Wave;

%% Discretization
%method = QDG;
%method = DGSEM;
%method = FR('DG');
%method = FR('Ga');
%method = FR('LumpLo');
%method = FR(1e-2);
method = DGIGA(4);
%method = DGIGA_AFC(4);

%% Initial condition
FUN = @(x) combinedIC(x);

%% Grid
mesh = Mesh(linspace(L(1),L(2),Ne+1),p,method);

%% Limiter
limiter = [];
%limiter = Limiter.TVBM(eqn,0);
%limiter = Limiter.TVBM(eqn,20);
%limiter = Limiter.Biswas(eqn);
%limiter = Limiter.Burbeau(eqn);
%limiter = Limiter.Krivodonova(eqn); % best in this case
%limiter = Limiter.Wang(eqn);
if isa(method,'DGIGA_AFC')
    limiter = Limiter.AFC(eqn);
end

%% Initial condition projection
% method.interpolate(mesh,limiter,FUN);
method.project(mesh,limiter,FUN);

%% Time-integration
timeIntegrator = SSP_RK3(0,tEnd,eqn,limiter,CFL,dt);
norms0 = [mesh.getSolutionNorm(1:2),mesh.getTotalVariation,mesh.getErrorNorm(FUN,1:2)];
tic
timeIntegrator.launch(mesh,iterSkip,@(t,x) FUN(x));
fprintf(1,'...done. (%g s)\n',toc)

%% Postprocessing
norms = [mesh.getSolutionNorm(1:2),mesh.getTotalVariation,mesh.getErrorNorm(FUN,1:2)];
rows = {'Solution (L2)','Solution (TV)','Error (L2)'};
cols = {'Norm','Start','End'};
eqn.displayData(rows,cols,norms0,norms)

%% Auxiliary functions
function y = constantIC(x)
% Initial condition that models a smooth displacement pulse.
y = zeros(2,length(x));
y(1,:) = 1;
end

function y = pickIC(x)
% Initial condition that models a smooth displacement pulse.
y = zeros(2,length(x));
y(1,:) = exp(-32*(x).^2);
end

function y = hammerIC(x)
% Initial condition that models a sudden perturbation in displacement rate 
% on a portion of an initially still medium.
y = zeros(2,length(x));
y(2,:) = heaviside(x-2/3+1) - heaviside(x-2*2/3+1);
end

function y = standingIC(x)
% Initial condition that will generate a smooth standing wave.
y = zeros(2,length(x));
y(2,:) = sin(2*x*pi);
end

function y = combinedIC(x)
% Initial condition that models the combination of a 'pick' (smooth pulse
% in displacement) and a 'hammer strike' (sudden change in displacement 
% rate) over non-overlapping regions of the domain.
y = zeros(2,length(x));
y(1,:) = Functions.gauss(x,0,.5);
y(2,:) = Functions.jump(x,.5,1);
end