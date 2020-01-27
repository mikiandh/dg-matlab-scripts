close all
clear
clc

%% Description
% Plots stability regions of various ODE time-integration schemes
%
%% Grid
% Specify x range and number of points
x0 = -6;
x1 =  2;
Nx = 401;
% Specify y range and number of points
y0 = -6;
y1 =  6;
Ny = 601;
% Construct mesh
xv    = linspace(x0,x1,Nx);
yv    = linspace(y0,y1,Ny);
[x,y] = meshgrid(xv,yv);
% Calculate z
z = x + 1i*y;

%% Time-scheme growth factor magnitudes
EF = @(z) abs(1 + z); % 1st order Euler
RK2 = @(z) abs(1 + z + 1/2*z.^2); % 2nd order Runge-Kutta
RK4 = @(z) abs(1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4); % 4th order Runge-Kutta
SSP_RK3_3 = @(z) abs(1 + z + 1/2*z.^2 + 1/6*z.^3); % 3rd order, 3 stage, SSP Runge-Kutta
SSP_RK4_10 = @(z) abs(z.^10/251942400 + z.^9/4199040 + z.^8/155520 + z.^7/9720 + (7*z.^6)/6480 + (17*z.^5)/2160 + z.^4/24 + z.^3/6 + z.^2/2 + z + 1);

SCHEMES = {EF RK2 RK4 SSP_RK3_3 SSP_RK4_10};
NAMES = {'Euler' 'RK2' 'RK4' 'SSP RK3(3)' 'SSP RK4(10)'};
COLORS = {'b','r','g','c','m','y'};

%% Stability contours
figure(1)
hold on
for i = 1:numel(SCHEMES)
contour(x,y,SCHEMES{i}(z),[1 1],COLORS{i},'DisplayName',NAMES{i});
end
hold off
axis([x0,x1,y0,y1]);
axis('equal');
xlabel('Real \lambda\Delta t');
ylabel('Imaginary \lambda\Delta t');
title('|g(\lambda\Delta t)|')
grid on;
legend('-DynamicLegend','Location','Best');