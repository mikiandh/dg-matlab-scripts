clc
clear

% This script test the exact solution to Burgers' equation with "sawtooth"
% initial condition.

%% Setup
a = 0;
t = 0:.1:2; % sample instants
x = linspace(-1,1,600); % solution sample points
nu = 1e-4; % artificial viscosity (Cole-Hopf solver)

%% Solutions
u = {
    sawtoothBurgersExact(t,x,a) 'Analytical'
    Burgers.MOC_robust(t,x,@(x) a+1-abs(x),[-1 1]) 'M.O.C.'
    burgersColeHopf(t,x,@(x) a+1-abs(x),nu,[-1 1]) sprintf('Cole-Hopf, %g',nu)
    };

%% Plot
for j = size(u,1):-1:1
    h(j) = plot(x,u{j,1}(1,:),'DisplayName',u{j,2});
    hold on
end
hold off
xlabel('x')
ylabel('u')
title(sprintf('t = %g',t(1)))
legend(h(end:-1:1))
axis manual
for i = 1:numel(t)
    for j = 1:numel(h)
        h(j).YData = u{j,1}(i,:);
    end
    title(sprintf('t = %g',t(i)))
    if i < numel(t) && ~mod(t(i),1)
        pause
    else
        pause()
    end
end