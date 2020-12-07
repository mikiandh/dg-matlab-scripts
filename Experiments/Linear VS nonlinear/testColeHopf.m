% This script test the Cole-Hopf transform solver for Burgers equation.

%% Setup
nu = 1e-3; % viscosity
t = 0:.05:.75; % sample instants
x = linspace(-1,1,500); % sample points
%%% f = @(x) 1+.1*sin(8*pi*x); % initial condition
f = @(x) exp(-9*pi/4*x.^2); % initial condition

%% Plot characteristic curves
x0 = linspace(x(1),x(end),50);
xi = f(x0).'.*t + x0.'; % trajectories (row: origin point; column: instant)
figure
yyaxis left
plot(xi,t,'-')
ylabel('t')
yyaxis right
plot(x,f(x),'-');
ylabel('f')
xlabel('x')

%% Plot initial condition
figure
plot(x,f(x),':k','DisplayName',sprintf('t = %g',0))
drawnow

%% Compute solution in time
y = burgersColeHopf(t,x,f,nu,x([1 end]));
z = Burgers.MOC_robust(t,x,f,x([1 end]));

%% Plot solution in time
hold on
h1 = plot(x,y(1,:),'DisplayName',sprintf('t = %g',t(1)));
h2 = plot(x,z(1,:),'DisplayName',sprintf('t = %g, inviscid',t(1)));
hold off
xlabel('x')
ylabel('u')
title(sprintf('Burgers (viscosity: %g); f = %s',nu,func2str(f)))
legend('Location','Best')
for i = 2:numel(t)
    h1.YData = y(i,:);
    h1.DisplayName = sprintf('t = %g, viscous',t(i));
    h2.YData = z(i,:);
    h2.DisplayName = sprintf('t = %g, inviscid',t(i));
    pause
end