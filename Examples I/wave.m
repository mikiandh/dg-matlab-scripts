clc
clear
%close all
%path(pathdef)

% This script solves the wave equation in 1D.

%% Dependencies
addpath('../Extra')
addpath('../Discretization')
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
Ne = 32; % number of elements (!-> for 1 element, use FIXED time-step size)
p = 2; % degree of the approximation space (per element)
L = [0 1]; % domain edges
tEnd = 2.0; % final simulation time
dt = [];
CFL = .1; % Courant number
iterSkip = 50;
    
%% Physics
eqn = Wave;

%% Discretization
%method = DG;
%method = DGSEM;
%method = FR('DG');
%method = FR('Ga');
%method = FR('LumpLo');
%method = FR(1e-2);
%method = DGIGA(28);
%method = DGIGA_AFC(30);
method = DG;

%% Initial condition
FUN = @(x) combinedIC(x);

%% Limiter
limiter = TVB;

%% Grid
mesh = Mesh(linspace(L(1),L(2),Ne+1),p,method,eqn);

%% Initial condition projection
% method.interpolate(mesh,limiter,FUN);
method.project(mesh,limiter,FUN);
% method.projectLumped(mesh,limiter,FUN);

%% Time-integration
solver = SSP_RK3(0,tEnd,CFL,dt,limiter);
norms0 = [mesh.getSolutionMass(1:2),mesh.getErrorNorm(FUN,2,1:2)];
tic
solver.launch(mesh,iterSkip,@(t,x) FUN(x));
fprintf(1,'...done. (%g s)\n',toc)

%% Postprocessing
norms = [mesh.getSolutionMass(1:2),mesh.getErrorNorm(FUN,2,1:2)];
rows = {'Solution (Mass)','Error (L2)'};
cols = {'Norm','Start','End'};
eqn.displayData(rows,cols,norms0,norms)

%% Auxiliary functions
function y = constantIC(x)
y = zeros(2,length(x));
y(2,:) = 1;
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