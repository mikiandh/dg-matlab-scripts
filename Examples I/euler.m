%clc
clear
%close all
%path(pathdef)

% This script solves the Euler equations.

%% Dependencies
addpath('../Extra')
addpath('../Discretization')
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Parameters
Ne = 300; % number of elements
p = 3; % degree of the approximation space (per element)
L = [-5 5]; % domain edges
tEnd = 1.8; % final simulation time
dt = [];
CFL = .1; % Courant number
iterSkip = 100;

%% Physics
eqn = Euler('transmissive');

%% Discretization
method = DG;

%% Limiter
limiter = BSB;

%% Grid
mesh = Mesh(linspace(L(1),L(2),Ne+1),p,method,eqn);
%mesh = Mesh(cosspace(L(1),L(2),Ne+1,2),p,method);
%mesh = Mesh([L(1) .05 .95 L(2)],p,method);

%% Initial condition/exact solution
FUN = @shuOsher;

%% Initial condition projection
FUN0 = @(x) FUN(0,x);
method.project(mesh,limiter,FUN0);

%% Time-integration
solver = SSP_RK3(0,tEnd,CFL,dt,limiter);
% norms0 = mesh.getSolutionMass(1:3);
norms0 = [mesh.getSolutionMass(1:3),mesh.getErrorNorm(FUN0,2,1:3)];
tic
solver.launch(mesh,iterSkip,FUN);
disp(['   ...done. (' num2str(toc) ' s)'])

%% Postprocessing
% norms = mesh.getSolutionMass(1:3);
% rows = {'Solution (mass)'};
norms = [mesh.getSolutionMass(1:3),mesh.getErrorNorm(FUN0,2,1:3)];
rows = {'Solution (mass)','Error (L2)'};
cols = {'Norm','Start','End','Increase'};
eqn.displayData(rows,cols,norms0,norms,norms-norms0)

%% Exact solutions and/or initial conditions
% Simple density jump:
function y = jump(t,x)
[r,u,p] = riemannEulerExact(t,x,2,0,1,1,0,1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Sod's shock tube:
function y = sod(t,x)
[r,u,p] = riemannEulerExact(t,x,1,0,1,0.125,0,0.1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Leveque 2002 (designed to avoid non-positive density due to undershoots):
function y = leveque(t,x)
[r,u,p] = riemannEulerExact(t,x,3,0.9,3,1,0.9,1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Toro (2009)'s collection (table 11.1):
function y = toro1(t,x)
% tEnd = 0.2
[r,u,p] = riemannEulerExact(t,x,1,0.75,1,0.125,0,0.1,0.3);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro2(t,x)
% tEnd = 0.15
[r,u,p] = riemannEulerExact(t,x,1,-2,0.4,1,2,0.4,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro2bis(t,x)
% Adapted to avoid non-positive density due to undershoots.
% tEnd = 0.15
[r,u,p] = riemannEulerExact(t,x,1,-.5,0.4,1,.5,0.4,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro3(t,x)
% tEnd = 0.012
[r,u,p] = riemannEulerExact(t,x,1,0,100,1,0,0.01,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro4(t,x)
% tEnd = 0.035
[r,u,p] = riemannEulerExact(t,x,5.99924,19.5975,460.894,5.99242,-6.19633,46.095,0.4);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro5(t,x)
% tEnd = 0.012
[r,u,p] = riemannEulerExact(t,x,1,-19.59745,1000,1,-19.59745,0.01,0.8);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Woodward and Colella:
function y = woodwardColella(~,x)
% tEnd = 0.038
y = ones(3,length(x));
y1 = Euler.primitiveToState([1 0 1000]');
y2 = Euler.primitiveToState([1 0 0.01]');
y3 = Euler.primitiveToState([1 0 100]');
y(:,x < 0.1) = y1.*y(:,x < 0.1);
y(:,x >= 0.1 & x < 0.9) = y2.*y(:,x >= 0.1 & x < 0.9);
y(:,x > 0.9) = y3.*y(:,x > 0.9);
end
% Shu and Osher (mimicks shock-turbulence interaction):
function y = shuOsher(~,x)
% tEnd = 1.8
% L = [-5 5]
y = ones(3,length(x));
for i = 1:length(x)
    if x(i) < -4
        y(:,i) = Euler.primitiveToState([3.857143 2.629369 10.33333]');
    else
        y(:,i) = Euler.primitiveToState([1+0.2*sin(5*x(i)) 0 1]');
    end
end
end
% Lax:
function y = lax(t,x)
% tEnd = 1.3
% L = [-5 5]
[r,u,p] = riemannEulerExact(t,x,0.445,0.698,3.528,0.5,0,0.571,0);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Smooth sine wave:
function y = sineSpeed(~,x)
y = ones(3,length(x));
y(2,:) = y(2,:).*sin(2*pi*x);
end