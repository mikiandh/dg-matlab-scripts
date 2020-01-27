clc
clear
%close all
%path(pathdef)

% This script performs a numerical modified wavenumber analysis for a given
% discretization scheme. Procedure adapted from van den Abeele 2009 (PhD
% thesis), section B1.1. All eigenmodes are considered.
%
% The generating pattern considered is a

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')

%% Parameters
% Discretization, e.g. DGSEM(2), DGIGA(2,1), DGIGA([-1 -1 0 1 1]):
disc = DGIGA(1,1,inf);
% Upwinding ratio for Riemann fluxes, from -1 (downwind) to 1 (upwind):
upwind = 1;
% Number of patches (only affects resolution):
Nx = 60;

%% Modified wavenumber analysis
[kMod,k] = MWA_eigen_full(Nx,disc.degree,disc,upwind);

%% Modified wavenumber plots
MWA_eigen_full_GUI(disc,k,kMod)