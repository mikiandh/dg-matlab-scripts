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

%% Example 2, Cockburn & Shu 1989 (III)
%
%  Euler's equation with transmissive BCs. Initial condition is Lax's
%  shock tube.
%
% Set up:
eqn = Euler('transmissive','HLLC');
limiter = Limiter.TVBM(eqn,0);
space = DG;
time = SSP_RK3(0,1.3,eqn,limiter,.1,[]);
fun = @lax;
fun0 = @(x) fun(0,x);
iterSkip = 100;
% Get exact solution:
x = linspace(-5,5,1000);
aux = fun(1.3,x);
y(3,:) = aux(1,:);
z{3} = sprintf('Exact');
% Get computed solution:
i = 0;
Ne = 100;
for p = [1 2]
    i = i+1;
    mesh = Mesh(linspace(-5,5,Ne+1),p,space);
    space.project(mesh,limiter,fun0);
    time.timeNow = 0;
    time.launch(mesh,iterSkip,fun);
    aux = mesh.sample(x);
    y(i,:) = aux(1,:);
    z{i} = sprintf('p = %d, Ne = %d',p,Ne);
end
% Plot solutions:
figure()
plot(x,y,'-')
xlabel('Position')
ylabel('Density')
setFancyPlotSettings1
title('DG, SSP-RK3; TVB, M2 = 300')
h = get(gca,'Children');
h(1).Color = 'k';
h(1).LineStyle = '-.';
legend(z{:},'Location','Best')

%% Example 3, Cockburn & Shu 1989 (III)
%
%  Euler's equation with transmissive BCs. Initial condition is Shu-Osher's
%  simple model of a shock-turbulence interaction.
%
% Set up:
eqn = Euler('transmissive','LLF');
limiter = Limiter.TVBM(eqn,200);
space = DG;
time = SSP_RK3(0,1.8,eqn,limiter,.1,[]);
fun = @shuOsher;
fun0 = @(x) fun(0,x);
iterSkip = 500;
% Get reference solution:
load example3
solutions(4).x = data(:,1);
solutions(4).y = data(:,2);
solutions(4).name = 'Reference';
% Get computed solution:
p = [1 2 2];
Ne = [200 100 200];
for i = 1:length(solutions)-1
    mesh = Mesh(linspace(-5,5,1+Ne(i)),[0 ones(1,Ne(i)-1)*p(i)],space);
    space.project(mesh,limiter,fun0);
    time.timeNow = 0;
    time.launch(mesh,iterSkip,fun);
    x = mesh.getDofLocations;
    y = mesh.sample(mesh.getDofLocations);
    solutions(i).x = x;
    solutions(i).y = y(1,:);
    solutions(i).name = sprintf('p = %d, Ne = %d',p(i),Ne(i));
end
% Plot solutions:
figure()
hold on
for aux = solutions
    plot(aux.x,aux.y,'.-','DisplayName',aux.name)
end
hold off
xlabel('Position')
ylabel('Density')
setFancyPlotSettings1
title('DG, SSP-RK3; TVB, M2 = 300')
h = get(gca,'Children');
h(1).Color = 'k';
h(1).LineStyle = '-.';
h(1).Marker = 'none';
legend('Location','Best')

%% Functions
% Lax:
function y = lax(t,x)
% tEnd = 1.3
% L = [-5 5]
[r,u,p] = riemannEulerExact(t,x,0.445,0.698,3.528,0.5,0,0.571,0);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
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