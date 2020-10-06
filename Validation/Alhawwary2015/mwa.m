clc
clear
close all
%path(pathdef)

% This script attempts to validate the energy-based definition of
% eigenmodes against Alhawwary and Wang, 2015.

%% Dependencies
addpath('../../Basis')

%% Central flux, dominant eigenmode
%
%  Reproduces (partially) figure 3b of Alhawwary and Wang (2015)  
%
load alhawwary3b.mat
figure(1)
subplot(2,1,1)
hold on
plot(DGp2beta0disp(:,1),DGp2beta0disp(:,2),'+k','DisplayName','Alhawwary and Wang (2015)')
hold off
DGSEM(2).displayDispDiss('beta',0,'order','energy')