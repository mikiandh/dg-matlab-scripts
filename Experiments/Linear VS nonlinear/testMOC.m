%clc
clear
close all

% This script test the M.O.C. solver for Burgers equation.

%% Setup
dt = .1; % sample frequency (in time)
x = linspace(-1,1,500); % solution sample points
x0 = x(1:10:end); % characteristic origin points

%% Test cases
f = {
    @(x) 1 + cos(pi*x/2).^2 % single solution, positive slopes
    @(x) -sin(pi*x) % single solution, opposing slopes
    @(x) -x % single solution, with multiplicity
    @(x) x % invalid solution (diverging lines)
    @(x) 1 + 0.*x % no solution (parallel lines)
    @(x) sin(pi*x) % multiple solutions
    @(x) heaviside(x-.25) % expansion fan
    @(x) 1-heaviside(x-.25) % shockwave
    @(x) 1.5-abs(x) % gentle breakdown
};

%% Solutions
for i = 1:numel(f)
    t = 0:dt:min(10,Burgers.getBreakingPoint(f{i},x([1 end])));
    subplot(1,2,1)
    yyaxis left
    h1(:,3) = plot(f{i}(x0).'.*t + x0.',t,'-'); % some characteristic lines
%     hold on
%     h1(:,2) = plot(vertcat(h1(:,3).XData)' + 2,vertcat(h1(:,3).YData)','--');
%     h1(:,1) = plot(vertcat(h1(:,3).XData)' - 2,vertcat(h1(:,3).YData)',':');
%     hold off
    ylabel('t')
    yyaxis right
    plot(x,f{i}(x),'-') % initial condition
    ylabel('f')
    xlabel('x')
    title(func2str(f{i}))
    subplot(1,2,2)
    legLabels = compose('%s, t = %g',["Robust" "Fast"]',t);
    tic
    u = Burgers.MOC_robust(t,x,f{i},x([1 end]));
    toc
    h3 = plot(x,u,':');
    [h3.DisplayName] = legLabels{1,:};
    tic
    for j = 1:numel(t)
        u(j,:) = Burgers.MOC(t(j),x,f{i},x([1 end]));
    end
    toc
    hold on
    h4 = plot(x,u,'-');
    hold off
    [h4.DisplayName] = legLabels{2,:};
    legend([h3; h4])
    if i < numel(f)
        pause
        clc
    end
end