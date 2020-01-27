clc, clear, %close all

%% Description
% Script that uses symbolic math to perform von Neumann linear stability 
% analysis for the advection equation, in 1D:
%
% \Delta t * (d/dt U_i) = -(\Delta t)/(\Delta x) * R_i
%
% Reference: Hirsch 2007

%% Symbolic variables
syms dx dt % cell and time-step sizes
syms U_k % amplitude of the k-th wavemode 
syms theta % phase angle (theta) = wavenumber (p_k) * cell size (dx)
syms id % cell id

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = pi; % domain size
N = 128; % number of cells
CFL = 0.8; % Courant number
lambda = 1; % advection speed
T = array2table([N L L/N pi/L pi*N/L CFL*L/N/lambda lambda CFL],'VariableNames',{'N','L','dx','k_min','k_max', 'dt', 'lambda', 'CFL'}) %#ok<NOPTS>

% Flux functions (flux on face id-1/2)
% >> Classic schemes (when combined with forward Euler) <<
LW = @(u,id) 0.5*lambda*(u(id-1) + u(id)) - 0.5*dt/dx*lambda^2*(u(id) - u(id-1)); % Lax-Wendroff
BW = @(u,id) lambda*u(id-1) + 0.5*lambda*(1-dt/dx*lambda)*(u(id-1) - u(id-2)); % Beam-Warming
% >> 1st order FD <<
U1 = @(u,id) lambda*u(id-1); % upwind
D1 = @(u,id) lambda*u(id); % downwind
% >> 2nd order FD <<
C2 = @(u,id) 0.5*lambda*(u(id-1) + u(id)); % centered
U2 = @(u,id) 0.5*lambda*(3*u(id-1) - u(id-2)); % upwind
% >> 4th order FD <<
C4 = @(u,id) -1/12*lambda*(u(id-2) - 7*u(id-1) - 7*u(id) + u(id+1)); % centered
%%%%%%%%%%%%
FLUX = C2;
%%%%%%%%%%%%

% Explicit ODE integrator:
TIME = ESSPRK(10,4); % SSP RK (stages, order)

% Try to plot g(z) isocontours?
flag_isocontours = 0; % set to FALSE if symbolic math gets stuck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution (single arbitrary wavemode)
syms u(id) % solution at cell id
u(id) = U_k*exp(1i*theta*id);

%% Assumptions
assume(dx,'positive')
assume(dt,'positive')
assume(U_k,'positive')
assume(theta,'real')
assume(theta >= 0 )
assume(id,'integer')
assume(id >= 0 )

%% Spatial operator
% Equation residual (i.e. net flux out of cell 'id'):
R = @(u) FLUX(u,id) - FLUX(u,id+1);
% Fourier symbol of the space operator:
z = simplify(R(u)/u);
% Space operator:
SPACE = @(u) 1/dx*R(u);

%% Temporal operator
% Future time-step solution:
uNext(id) = simplify(TIME(dt,u,SPACE));

%% Amplification factor
g = simplify(uNext/u);
g = simplifyFraction(g,'expand',1);

if flag_isocontours
    % Some necessary tricks... (for the LHS subplot)
    theta_z(theta) = finverse(z,theta); % theta as a function of z (independent variable => z)
    g_z = simplify(subs(g,theta,theta_z)); % g as a function of z
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF THE SYMBOLIC MATH PART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visualization of results
% Amplification factor vs. phase angle (i.e. diffusion and phase errors):
theta_h = linspace(-pi,pi,N); % theta, spanning half the phase space (other half is symmetric)
g_h = matlabFunction(g,'Vars',[id dx dt theta]);
g_h = g_h(0,L/N,CFL*L/N/lambda,theta_h);
g_abs = abs(g_h);
g_phase = -angle(g_h); % sign?
g_phase_exact = theta_h*CFL;

% Amplification factor vs. spatial fourier symbol:
z_h = matlabFunction(z,'Vars',[id dx dt theta]);
z_h = z_h(0,L/N,CFL*L/N/lambda,theta_h); % locus of the complex plane that corresponds to the phase angle under current conditions (e.g. CFL)
z_real = real(z_h);
z_imag = imag(z_h);
if flag_isocontours
    % Complex plane grid:
    y0 = min(-3.5,min(z_imag)-2);
    y1 = max(3.5,max(z_imag)+2);
    x1 = 1;
    x0 = x1 - (y1 - y0);
    [x,y] = meshgrid(linspace(x0,x1,350),linspace(y0,y1,350));
    z_grid = x + 1i*y; % spatial fourier symbol as an independent variable
    g_grid = matlabFunction(g_z,'Vars',[id dx dt theta]); % from symbolic expression to anonymous function (speeds up evaluation)
    g_abs_grid = abs(g_grid(0,L/N,CFL*L/N/lambda,z_grid)); % evaluate function
end

% Plots:
figure('Renderer', 'painters', 'Position', [25 100 1800 500])

subplot(1,3,1)
plot(z_real,z_imag,'k','LineWidth',2)
hold on
if flag_isocontours
    contour(x,y,g_abs_grid,0:0.2:1,'ShowText','on')
end
hold off
xlabel('Real(z)')
ylabel('Imag(z)')
legend('|G|(\Phi)','|G|(z)','Location','NorthWest')
setFancyPlotSettings2
axis equal
title(['CFL = ' num2str(CFL)])

subplot(1,3,2)
plot(theta_h,g_abs)
hold on
hline = refline([0 1]);
hline.Color = 'k';
hline.LineStyle = '--';
hold off
axis([0 pi 0 inf])
xticks([0 0.25*pi 0.5*pi 0.75*pi pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('Phase angle, \Phi (radians)')
ylabel('Magnitude of amplification factor, |G| (-)')
setFancyPlotSettings2
title(['CFL = ' num2str(CFL)])

subplot(1,3,3)
plot(theta_h,g_phase)
hold on
plot(theta_h,g_phase_exact,'--k')
hold off
xlim([0 pi])
xticks([0 0.25*pi 0.5*pi 0.75*pi pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
ylim([0 inf])
yticks([0 0.25*pi 0.5*pi 0.75*pi pi 1.25*pi 1.5*pi 1.75*pi 2*pi])
yticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
xlabel('Phase angle, \Phi (radians)')
ylabel('Phase of amplification factor, \angle G (radians)')
setFancyPlotSettings2
title(['CFL = ' num2str(CFL)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF THE SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%