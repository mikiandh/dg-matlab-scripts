clc
clear
%close all

% This script solves the following nonlinear, scalar, real-valued
% optimization problem:
%
% Given FR of degree 'p' and an explicit time scheme, find the correction
% parameter 'eta' that minimizes the L2-norm of the error with the exact
% dispersion relation.

%% Dependencies
addpath('../../Solver')
addpath('../../Basis')

%% Setup
p = [2 3 4 5];% 7 10 14 19];
switch 1
    case 1
        objFun = @(eta,p) objFun_peakPosition(FR({'eta',eta},p),.5);
        filename = 'fr_dispersion_peak50.dat';
    case 2
        objFun = @(eta,p) objFun_proportionalDissipation(FR({'eta',eta},p),1,true);
        filename = 'fr_dispersion_prop100.dat';
end
time = SSP_RK3;

%% Preprocess
try
    tblIn = readtable(filename);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(p);
tbl = array2table(zeros(I,7),'VariableNames',{'p','eta','c','exitFlag','relCost','A_T','CFL'});
for i = 1:I
    try
        options = optimset('PlotFcn',{
            @(varargin) plotFun_dispDiss(p(i),varargin{:})
        });
        [eta,badness,flag] = fminbnd(@(x) objFun(x,p(i)),0,5,options);
        basis = FR({'eta',eta},p(i));
        badness = badness/objFun(0,p(i)) - 1; % relative to DG
        tbl{i,:} = [p(i) basis.eta basis.c flag badness basis.getOrder time.optimizeCFL(basis)];
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
    catch me
       warning("Run %d (p = %d) failed with error:\n '%s'",i,p(i),getReport(me))
    end
end

%% Postprocess
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
writetable(tblOut,filename,'Delimiter','tab')