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
BE = @(z) abs(1./(1 - z)); % 1st order backward Euler
FE = @(z) abs(1 + z); % 1st order forward Euler
SSP_RK2 = @(z) abs(1 + z + 1/2*z.^2); % 2nd order Runge-Kutta
SSP_RK3 = @(z) abs(1 + z + 1/2*z.^2 + 1/6*z.^3); % 3rd order, 3 stage, SSP Runge-Kutta
RK4_4 = @(z) abs(1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4); % 4th order Runge-Kutta
SSP_RK4_5 = @(z) abs((265463990952622608150189972344896021762287709730189019808505249429955773138217*z.^5)/59285549689505892056868344324448208820874232148807968788202283012051522375647232 + (617557809265686329792262629871895994724553962772224133823548839147548176662523*z.^4)/14821387422376473014217086081112052205218558037201992197050570753012880593911808 + (2470231237062745899416026889164753911888427586318335394030467413516256206203135*z.^3)/14821387422376473014217086081112052205218558037201992197050570753012880593911808 + (1852673427797059486576953433340741858439650075155231662413332087780857683018485*z.^2)/3705346855594118253554271520278013051304639509300498049262642688253220148477952 + (231584178474632494065353036093698953526959801287985692919653483466816886494293*z)/231584178474632390847141970017375815706539969331281128078915168015826259279872 + 205688069665150948943867075817191687448919604710488653564680541/205688069665150755269371147819668813122841983204197482918576128);
SSP_RK4_10 = @(z) abs(z.^10/251942400 + z.^9/4199040 + z.^8/155520 + z.^7/9720 + (7*z.^6)/6480 + (17*z.^5)/2160 + z.^4/24 + z.^3/6 + z.^2/2 + z + 1);

SCHEMES = {BE FE SSP_RK2 SSP_RK3 RK4_4 SSP_RK4_5 SSP_RK4_10};
NAMES = {'BE' 'FE' 'SSP RK2(2)' 'SSP RK3(3)' 'RK4(4)' 'SSP RK4(5)' 'SSP RK4(10)'};
COLORS = {'b','r','g','c','m','y','b'};

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