clc
clear
close all

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
p = logspacei(2,19,8);
objFun = @(eta,p) objFun_Asthana2015(FR({'eta',eta},p));
time = SSP_RK3;
export = struct('name','fr_dispersion',...
    'dat',true,...
    'fig',false,...
    'gif',false,... incompatible with parfor, sorry
    'tikz',false);

%% Preprocess
try
    tblIn = readtable([export.name '.dat']);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(p);
tbl = array2table(zeros(I,15),'VariableNames',{
    'p'
    'eta'
    'c'
    'Goodness'
    'relGoodness'
    'DispDiss'
    'relDispDiss'
    'Order'
    'relOrder'
    'Resol'
    'relResol'
    'Cutoff'
    'relCutoff'
    'CFLRK3'
    'relCFLRK3'
    });
parfor i = 1:I
    try
        options = optimset('PlotFcn',{
            @(varargin) plotFun_dispDiss(varargin{:},p(i),i,export)
        });
        [eta,badness] = fminbnd(@(x) objFun(x,p(i)),-1,5,options);
        optimum = FR({'eta',eta},p(i));
        baseline = DGSEM(p(i));
        [~,~,optimumDispDiss] = optimum.getDispDissRatio;
        [~,~,baselineDispDiss] = baseline.getDispDissRatio;
        tbl{i,:} = [
            p(i)
            optimum.eta
            optimum.c
            -badness
            -objFun_Asthana2015(baseline)
            optimumDispDiss
            baselineDispDiss
            optimum.getOrder
            baseline.getOrder
            optimum.getResolvingWavenumber/optimum.basisCount/pi
            baseline.getResolvingWavenumber/baseline.basisCount/pi
            optimum.getCutoffWavenumber/optimum.basisCount/pi
            baseline.getCutoffWavenumber/baseline.basisCount/pi
            time.optimizeCFL(optimum)
            time.optimizeCFL(baseline)
        ]'; %#ok<PFBNS>
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
        drawnow
        %%%
        if export.fig
            saveas(gcf,[export.name '_' num2str(i) '.fig']);
        end
        if export.tikz
            cleanfigure
            matlab2tikz([export.name '_' num2str(i) '.tikz'],'showInfo',false)
        end
        %%%
    catch me
       warning("Run %d (p = %d) failed with error:\n '%s'",i,p(i),getReport(me))
    end
end

%% Postprocess
if ~isempty(tbl)
    tbl{:,5:2:end} = (tbl{:,4:2:end-1} - tbl{:,5:2:end})./abs(tbl{:,5:2:end}); % relative changes over baseline
end
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[export.name '.dat'],'Delimiter','tab')
end