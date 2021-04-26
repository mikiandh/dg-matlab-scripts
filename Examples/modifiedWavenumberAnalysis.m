clc
clear

% This script showcases the main variations of the modified wavenumber
% analysis that this code supports.

%% Discretization
space = DGSEM(2);
time = SSP_RK3;

%% Dispersion and dissipation relations (Bloch waves)
figure
space.displayModifiedWavenumbers

%% Dispersion and dissipation relations (physical mode)
figure
space.displayDispDiss('upwind',1)
space.displayDispDiss('upwind',0)

%% Dispersion and dissipation relations (all modes, relative energy shown)
figure
space.displayDispDissEnergy

%% Relative energy of each eigenmode
figure
space.displayEnergy

%% Dispersion and dissipation combined-mode errors (phase and amplitude)
figure
space.displayAngAmp([1 10 100])

%% Scalar metrics
space.getOrder % theoretical order of accuracy
space.getCutoffWavenumber % dissipation cutoff
space.getResolvingWavenumber % dissipation + dispersion cutoff