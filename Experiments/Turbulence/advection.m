clc
clear
%close all

% This script computes a posteriori modified wavenumbers for the Advection
% equation, using Stefan et al. (2014)'s approach.

%% Input
basis = DGSEM(1);
N = 5; % number of wavemodes
K = 10; % number of patches
M = 10; % number of sample points (mesh-wide, preferably a power of 2)

%% Setup
kMod = zeros(basis.basisCount,N);
mesh = Mesh(basis,[-1 1],Periodic(2),K);
solver = SSP_RK3(Advection,[0 0]);

%% Wavemode loop
for n = 1:N
    dx = mesh.elements(1).dx;
    J = mesh.bases.basisCount;
    k = 2*pi/(K*dx);
    fun = @(x) 1 + .1*sin(k*n*x);
    solver.initialize(mesh,'initialCondition',fun)
    close % automatic monitor figure
    x = linspace(mesh.edges([1 end]).coord,M);
    y0 = fun(x);
    y = mesh.sample(x);
    z = mesh.sampleResidual(x);
    fx = [0:M/2-1 -M/2:-1];
    fy0 = fft(y0);
    fy = fft(y);
    fz = fft(z);
    subplot(1,2,1)
    plot(x,fun(x),x,y,x,z)
    xlabel('x')
    legend('q','q^h','\partial_t q^h','Location','northoutside')
    subplot(1,2,2)
    bar(fx,[abs(fy0); abs(fy); abs(fz)].'/M)
    hold on
    plot([n n],[0 1],'--k')
    hold off
    legend('F(q)','F(q^h)','F(\partial_t q^h)','\kappa_{exact}','Location','northoutside')
    xlabel('n')
    pause
end