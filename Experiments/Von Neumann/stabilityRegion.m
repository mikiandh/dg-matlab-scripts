close all
clear
clc

% This script shows the stability regions of various amplification factors,
% each in its own figure, and all |G| = 1 contours together.

%% Time schemes
SCHEMES = {...
    %@(z) 1./(1 - z); % 1st order backward Euler
	@(z) SSP_RK1.amplificationFactorFun(z) % 1st order forward Euler
	@(z) SSP_RK2.amplificationFactorFun(z) % 2nd order Runge-Kutta
	@(z) SSP_RK3.amplificationFactorFun(z) % 3rd order, 3 stage, SSP Runge-Kutta
	%@(z) 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4 % 4th order Runge-Kutta
	@(z) SSP_RK4_5.amplificationFactorFun(z)
	@(z) SSP_RK4_10.amplificationFactorFun(z)
};
NAMES = {...
    %'Backward Euler'
    'SSP-RK1(1)'
    'SSP-RK2(2)'
    'SSP-RK3(3)'
    %'RK4(4)'
    'SSP-RK4(5)'
    'SSP-RK4(10)'
};
COLORS = distinguishable_colors(numel(SCHEMES),{'w', 'k'});

%% Grid
[x,y] = meshgrid(linspace(-16,2,250),linspace(-8,8,250));
z = x + 1i*y; % sample points

%% Stability regions
for i = 1:numel(SCHEMES)
    figure(i)
    contourf(x,y,abs(SCHEMES{i}(z)),[0 1]);
    colormap([COLORS(i,:); 1 1 1])
    % Coordinate axes:
    hold on
    plot([x(1) x(end)],[0 0],'k--')
    plot([0 0],[y(1) y(end)],'k--')
    hold off
    axis('equal');
    xlabel('Real \lambda\Delta t');
    ylabel('Imaginary \lambda\Delta t');
    title(NAMES{i})
    grid on;
end

%% Stability contours
figure
hold on
for i = 1:numel(SCHEMES)
    contour(x,y,abs(SCHEMES{i}(z)),[1 1],'LineColor',COLORS(i,:),'DisplayName',NAMES{i});
end
% Coordinate axes:
hasbehavior(plot([x(1) x(end)],[0 0],'k--'),'legend',false);
hasbehavior(plot([0 0],[y(1) y(end)],'k--'),'legend',false);
hold off
axis('equal');
xlabel('Re(z)');
ylabel('Im(z)');
legend('-DynamicLegend','Location','Best');
grid on;