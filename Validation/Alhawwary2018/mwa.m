clc
clear
close all
%path(pathdef)

% This script attempts to validate stuff against Alhawwary and Wang, 2018.

%% Dependencies
addpath('../../Basis')

%% Central flux, dominant eigenmode
%
%  Reproduces (partially) figure 3b of Alhawwary and Wang (2018)  
%
load alhawwary3b.mat
figure(1)
subplot(2,1,1)
hold on
plot(DGp2beta0disp(:,1),DGp2beta0disp(:,2),'+k','DisplayName','Alhawwary and Wang (2018)')
hold off
DGSEM(2).displayDispDiss('upwind',0,'criterion','energy')

%% Combined-mode analysis 1
%
% Reproduces figure 4a of Alhawwary and Wang (2018)
%
load alhawwary4a.mat
load alhawwary4b.mat
figure(2)
subplot(2,1,1)
hold on
plot(phase_DGp2_t01_combined(:,1),10.^phase_DGp2_t01_combined(:,2),'-.k','DisplayName','Ref.; t = 0.1, combined')
plot(phase_DGp2_t01_dominant(:,1),10.^phase_DGp2_t01_dominant(:,2),'o-.k','DisplayName','Ref.; t = 0.1, dominant')
plot(phase_DGp2_t1_combined(:,1),10.^phase_DGp2_t1_combined(:,2),'-b','DisplayName','Ref.; t = 1, combined')
plot(phase_DGp2_t1_dominant(:,1),10.^phase_DGp2_t1_dominant(:,2),'-^b','DisplayName','Ref.; t = 1, dominant')
plot(phase_DGp2_t10_combined(:,1),10.^phase_DGp2_t10_combined(:,2),'--c','DisplayName','Ref.; t = 10, combined')
plot(phase_DGp2_t10_dominant(:,1),10.^phase_DGp2_t10_dominant(:,2),'s--c','DisplayName','Ref.; t = 10, dominant')
hold off
subplot(2,1,2)
hold on
plot(amplitude_DGp2_t01_combined(:,1),amplitude_DGp2_t01_combined(:,2),'-.k','DisplayName','Ref.; t = 0.1, combined')
plot(amplitude_DGp2_t01_dominant(:,1),amplitude_DGp2_t01_dominant(:,2),'o-.k','DisplayName','Ref.; t = 0.1, dominant')
plot(amplitude_DGp2_t1_combined(:,1),amplitude_DGp2_t1_combined(:,2),'-b','DisplayName','Ref.; t = 1, combined')
plot(amplitude_DGp2_t1_dominant(:,1),amplitude_DGp2_t1_dominant(:,2),'-^b','DisplayName','Ref.; t = 1, dominant')
plot(amplitude_DGp2_t10_combined(:,1),amplitude_DGp2_t10_combined(:,2),'--c','DisplayName','Ref.; t = 10, combined')
plot(amplitude_DGp2_t10_dominant(:,1),amplitude_DGp2_t10_dominant(:,2),'s--c','DisplayName','Ref.; t = 10, dominant')
hold off
DGSEM(2).displayDispDissErrors(0.1)
DGSEM(2).displayDispDissErrors(0.1,'mode','dominant')
DGSEM(2).displayDispDissErrors(1)
DGSEM(2).displayDispDissErrors(1,'mode','dominant')
DGSEM(2).displayDispDissErrors(10)
DGSEM(2).displayDispDissErrors(10,'mode','dominant')