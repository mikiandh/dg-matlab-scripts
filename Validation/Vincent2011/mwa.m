clc
clear
close all
%path(pathdef)

% This script attempts to validate FR/CPR against Vincent et al., 2011
% (both papers from that year).

%% Dependencies
addpath('../../Basis')

%% Left correction functions
%
%  Reproduces figure 1b of Vincent et al. (2011, JSC)
%
load vincentJSC.mat
figure(1)
hold on
x = linspace(-1,1);
plot(x,FR.getCorrectionFunctionAndDerivative(FR('min',3).c,3,x))
plot(x,FR.getCorrectionFunctionAndDerivative(0,3,x))
plot(x,FR.getCorrectionFunctionAndDerivative(FR('Ga',3).c,3,x))
plot(x,FR.getCorrectionFunctionAndDerivative(FR('LumpLo',3).c,3,x))
plot(x,FR.getCorrectionFunctionAndDerivative(FR('max',3).c,3,x))
scatter(cmin_2(:,1),cmin_2(:,2),'+k')
scatter(c0(:,1),c0(:,2),'+k')
scatter(cSD(:,1),cSD(:,2),'+k')
scatter(cHU(:,1),cHU(:,2),'+k')
scatter(cInf(:,1),cInf(:,2),'+k')
hold off

%% Dispersion, c < 0
%
%  Reproduces figure 1a of Vincent et al. (2011, JCP).
%
load vincent1a.mat
c_min = FR('min',3).c*2;
bases = {
    FR(3*c_min/4,3)
    FR(2*c_min/4,3)
    FR(1*c_min/4,3)
    FR(0*c_min/4,3)
};
figure(2)
hold on
for i = 1:numel(bases)
    [z,k] = bases{i}.getFourierFootprint;
    plot(k(k >= 0),-imag(z(1,k >= 0)))
end
plot(k(k >= 0),k(k >= 0),'k--')
scatter(c3_4(:,1),c3_4(:,2),'+k')
scatter(c2_4(:,1),c2_4(:,2),'+k')
scatter(c1_4(:,1),c1_4(:,2),'+k')
scatter(c0_4(:,1),c0_4(:,2),'+k')
hold off

%% Dispersion, c > 0
%
%  Reproduces figure 1b of Vincent et al. (2011, JCP).
%
load vincent1b.mat
bases = {
    FR('DG',3)
    FR('Ga',3)
    FR('LumpLo',3)
};
figure(3)
hold on
for i = 1:numel(bases)
    [z,k] = bases{i}.getFourierFootprint;
    plot(k(k >= 0),-imag(z(1,k >= 0)))
end
plot(k(k >= 0),k(k >= 0),'k--')
scatter(cSD(:,1),cSD(:,2),'+k')
scatter(cHU(:,1),cHU(:,2),'+k')
hold off

%% Dispersion, c >> 0
%
%  Reproduces figure 1c of Vincent et al. (2011, JCP).
%
load vincent1c.mat
cHU = FR('LumpLo',3).c;
bases = {
    FR(cHU,3)
    FR(2*cHU,3)
    FR('max',3)
};
figure(4)
hold on
for i = 1:numel(bases)
    [z,k] = bases{i}.getFourierFootprint;
    plot(k(k >= 0),-imag(z(1,k >= 0)))
end
plot(k(k >= 0),k(k >= 0),'k--')
scatter(c1HU(:,1),c1HU(:,2),'+k')
scatter(c2HU(:,1),c2HU(:,2),'+k')
scatter(cInf(:,1),cInf(:,2),'+k')
hold off