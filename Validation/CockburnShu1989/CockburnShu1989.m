clc
clear
close all
%path(pathdef)

% This script tests TVB-limited modal DG a la Cockburn & Shu, 1989.
% Validates results against reported ones for various equations, boundary 
% and initial conditions.

%% Dependencies
addpath('../../Extra')
addpath('../../Discretization')
addpath('../../Limiters')
addpath('../../Physics')
addpath('../../Solver')
addpath('../../Grid')
addpath('../../Math')

%% Example 1, Cockburn & Shu 1989 (II)
%
%  Burger's equation with periodic BCs. Initial condition is smooth, but
%  develops into a shock at t = 2/pi.
%
% Set up:
load example1
eqn = Burgers;
method = DG;
limiter = Limiter.TVBM(eqn,pi^2/3);
fun0 = @(x) .25 + .5*sin(pi*x);
fun = @(t,x) smoothBurgersExact(t,x,fun0);
dt = 1e-3;
p = 2;
iterSkip = 100;
% Ne = 20, p = 3 0 < t < 0.3:
tic
mesh = Mesh(linspace(-1,1,21),p,method);
method.project(mesh,limiter,fun0);
timeIntegrator = SSP_RK3(0,.3,eqn,limiter,[],dt);
timeIntegrator.launch(mesh,iterSkip,fun);
toc
errors = table(2/mesh.elementCount,...
    mesh.getErrorNorm(@(x) fun(.3,x),1),...
    nan,...
    mesh.getErrorNorm(@(x) fun(.3,x),2),...
    nan);
% 0.3 < t < 2/pi:
timeIntegrator.timeStop = 2/pi;
timeIntegrator.launch(mesh,iterSkip,fun);
toc
% pi/2 < t < 1.1:
timeIntegrator.timeStop = 1.1;
timeIntegrator.launch(mesh,iterSkip,@(t,x) nan);
toc
% Ne = 40, p = 3, 0 < t < 0.3:
tic
mesh = Mesh(linspace(-1,1,41),p,method);
method.project(mesh,limiter,fun0);
timeIntegrator = SSP_RK3(0,.3,eqn,limiter,[],dt);
timeIntegrator.launch(mesh,iterSkip,fun);
toc
errors = [errors; table(2/mesh.elementCount,...
    mesh.getErrorNorm(@(x) fun(.3,x),1),...
    nan,...
    mesh.getErrorNorm(@(x) fun(.3,x),2),...
    nan)];
% 0.3 < t < 2/pi:
timeIntegrator.timeStop = 2/pi;
timeIntegrator.launch(mesh,iterSkip,fun);
toc
hold on
plot(data1(:,1),data1(:,2),'sk')
hold off
% pi/2 < t < 1.1:
timeIntegrator.timeStop = 1.1;
timeIntegrator.launch(mesh,iterSkip,@(t,x) nan);
toc
hold on
plot(data3(:,1),data3(:,2),'sk')
hold off
% Ne = 80, p = 3, 0 < t < 0.3:
tic
mesh = Mesh(linspace(-1,1,81),p,method);
method.project(mesh,limiter,fun0);
timeIntegrator = SSP_RK3(0,.3,eqn,limiter,[],dt);
timeIntegrator.launch(mesh,iterSkip,fun);
toc
errors = [errors; table(2/mesh.elementCount,...
    mesh.getErrorNorm(@(x) fun(.3,x),1),...
    nan,...
    mesh.getErrorNorm(@(x) fun(.3,x),2),...
    nan)];
% 0.3 < t < 2/pi:
timeIntegrator.timeStop = 2/pi;
timeIntegrator.launch(mesh,iterSkip,fun);
toc
hold on
plot(data2(:,1),data2(:,2),'sk')
hold off
% pi/2 < t < 1.1:
timeIntegrator.timeStop = 1.1;
timeIntegrator.launch(mesh,iterSkip,@(t,x) nan);
toc
hold on
plot(data4(:,1),data4(:,2),'sk')
hold off
% Order of accuracy (smooth):
errors.Properties.VariableNames = {'dx', 'ErrorL1', 'OrderL1', 'ErrorL2', 'OrderL2'};
getOrder = @(h,e) diff(log(e))./diff(log(h));
errors.OrderL1(2:end) = getOrder(errors.dx,errors.ErrorL1);
errors.OrderL2(2:end) = getOrder(errors.dx,errors.ErrorL2);
disp(errors)
