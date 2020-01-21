clc, clear, close all

%% Description
% Script that tests the exact Riemann solver for the 1D Euler equations

%% Intial conditions
% Interface:
x0 = 0.5;
% Left state:
rL = 1;
uL = 2;
pL = 1;
% Right state:
rR = 1;
uR = -2;
pR = 1;

%% Space
x = linspace(0,1,1000);

%% Time
t = linspace(0,0.1,10);

%% Solution
[r,u,p] = riemannEulerExact(t,x,rL,uL,pL,rR,uR,pR,x0);

%% Plot in time
figure('Renderer', 'painters', 'Position', [25 100 900 800])
for i = 1:length(t)
    subplot(3,1,1)
    plot(x,r(i,:))
    title(['t = ' num2str(t(i))])
    xlabel('x')
    ylabel('\rho')
    setFancyPlotSettings3
    subplot(3,1,2)
    plot(x,u(i,:))
    xlabel('x')
    ylabel('u')
    setFancyPlotSettings3
    subplot(3,1,3)
    plot(x,p(i,:))
    xlabel('x')
    ylabel('p')
    setFancyPlotSettings3
    pause(0.125)
end