clc
clear
%close all

% This script solves the following nonlinear, scalar, real-valued
% optimization problem:
%
% Given FR of degree 'p' and an explicit time scheme, find the correction
% parameter 'eta' that maximizes the maximum stable Courant number,
% 'CFL_max', of the combination.

%% Dependencies
addpath('../../Solver')
addpath('../../Basis')

%% Setup
p = 6:20;
time = SSP_RK3;
objFun = @(eta,p) -time.optimizeCFL(FR({'eta',eta},p));
filename = 'fr_stability.dat';

%% Preprocess
try
    tblIn = readtable(filename);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end
tbl = array2table(zeros(numel(p),5),'VariableNames',{'p','eta','c','CFL','exitFlag'});

%% Minimization
I = numel(p);
parfor i = 1:I
    try
        [eta,CFL,flag] = fminbnd(@(x) objFun(x,p(i)),0,10);
        basis = FR({'eta',eta},p(i));
        tbl{i,:} = [p(i) basis.eta basis.c -CFL flag];
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
    catch
        warning('Run %d (p = %d) failed.',i,p(i))
    end
end

%% Postprocess
tblOut = sortrows([tblIn; tbl],'p');
clc, disp(tblOut)
writetable(tblOut,filename,'Delimiter','tab')