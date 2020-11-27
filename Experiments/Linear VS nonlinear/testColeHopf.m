% This script test the Cole-Hopf transform solver for Burgers equation.

%% Setup
nu = 1e-3; % viscosity
t = 0:.05:4; % sample instants
x = linspace(-1,1,500); % sample points
%%% f = @(x) 1+.1*sin(2*pi*x); % initial condition
f = @(x) exp(-9*pi/4*x.^2); % initial condition

%% Plot initial condition
plot(x,f(x),':k','DisplayName',sprintf('t = %g',0))
drawnow

%% Compute solution in time
y = burgersColeHopf(t,x,f,nu,x([1 end]));

%% Plot solution in time
hold on
h = plot(x,y(1,:),'DisplayName',sprintf('t = %g',t(1)));
hold off
xlabel('x')
ylabel('u')
title(sprintf('Viscosity: %g',nu))
legend('Location','Best')
for i = 2:numel(t)
    h.YData = y(i,:);
    h.DisplayName = sprintf('t = %g',t(i));
    pause(.1)
end