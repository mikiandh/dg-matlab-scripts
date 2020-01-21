clc
clear
%close all
%path(pathdef)

% This script approximates the maximum time-step size bounds for a given 
% time and space semi-discretization. Procedure adapted from van den Abeele 
% 2009 (PhD thesis), section B3.
%
% The stability limit occurs when one point in the scaled Fourier footprint
% of the spatial scheme reaches the edge of the time scheme's stability 
% region - i.e. the amplification factor of the full (space and time)
% discretization first reaches unity as the time-step size decreases from 
% infinity.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Solver')

%% Temporal discretization
time = SSP_RK1;
% Temporal amplification factor function:
% G = @(z) 1./(1 - z); % implict RK1
G = time.amplificationFactorFun; % explicit SSP RK
% Courant number seed:
CFL = 10;

%% Spatial discretization
% Assumed uniform, on a periodic domain e.g. DGSEM(2), DGIGA(2,1),
% DGIGA([-1 -1 0 1 1]).
disc = DGIGA(1,0,inf);
% Upwinding ratio for Riemann fluxes, from -1 (downwind) to 1 (upwind):
upwind = 1;
% Number of patches (only affects resolution):
Nx = 15;

%% Modified wavenumber analysis
[kMod,k] = MWA_eigen_full(Nx,disc.degree,disc,upwind); % kMod = 1i*eigenvalue
% kMod = kMod(1,:); % consider only some eigenmodes
% kMod = kMod(:,end/4:end/4:end); % consider only some wavenumbers

%% Maximum stable Courant number search
theta = -1i*kMod; % eigenvalues
% Nonlinear constrained optimization:
problem.options = optimoptions('fmincon','Display','notify-detailed');
problem.solver = 'fmincon';
problem.objective = @(CFL) - sum(sum(abs(G(theta*CFL)),2)); % maximize some norm of amplitude factors, all eigenmodes at once
problem.x0 = CFL;
problem.lb = 0;
problem.nonlcon = @(CFL) nonlcon(abs(G(theta*CFL)));
tic
CFL = fmincon(problem);
toc

%% Verification plots
figure('Renderer', 'painters', 'Position', [100 100 1200 600])
k = horzcat(-flip(k),k); theta = horzcat(flip((theta').',2),theta); % recover negative wavenumbers (aesthetic reasons only)
% Full amplification factor (space + time) vs. wavenumber:
subplot(1,2,1)
g = abs(G(theta*CFL));
plot(k/disc.basisCount,g','-')
setFancyPlotSettings3
xlabel('$$\kappa \frac{\Delta x}{p + N_\sigma}$$','Interpreter','latex')
ylabel('$$|G(\varsigma\,\tilde{\Theta})|$$','Interpreter','latex')
% Fourier footprint (space) vs. stability region (time):
subplot(1,2,2)
[x,y] = meshgrid(linspace(-16,2,250),linspace(-8,8,250));
z = abs(G(x+1i*y));
contourf(x,y,z,[0 1]);
colormap([.8 1 1; 1 1 1]);
hold all
plot([x(1) x(end)],[0 0],'k--')
plot([0 0],[y(1) y(end)],'k--')
z = theta*CFL;
x = real(z);
y = imag(z);
[x,y] = points2contour(x(:),y(:),1,'cw');
x = [x x(1)];
y = [y y(1)];
plot(x,y,'-b') % join all modes in a single smooth contour
% plot(x',y','.') % plot each mode separately
hold off
setFancyPlotSettings3
% axis tight
axis equal
xlabel('$$\Re(\tilde{z})$$','Interpreter','latex')
ylabel('$$\Im(\tilde{z})$$','Interpreter','latex')
% Global title:
if disc.isHybrid
    Ns = disc.nonzeroSpanCount;
else
    Ns = 1;
end
title(sprintf('%s, p = %d, Ns = %d, Nx = %d \n CFL = %g \n',class(disc),disc.degree,Ns,Nx,CFL))

%% Constraint function
function [c,ceq] = nonlcon(g)
    c = max(max(g)) - 1;
    ceq = [];
end