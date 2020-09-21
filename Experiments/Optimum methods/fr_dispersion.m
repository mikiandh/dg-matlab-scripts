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
p = 5; % [2 3 4 5 7 10 14 19];
switch 0
    case 0
        objFun = @(eta,p) objFun_exactBeforeCutoff(FR({'eta',eta},p));
        name = 'fr_dispersion_exact';
    case 1
        objFun = @(eta,p) objFun_peakPosition(FR({'eta',eta},p),.5);
        name = 'fr_dispersion_peak50';
    case 2
        objFun = @(eta,p) objFun_proportionalDissipation(FR({'eta',eta},p),1,true);
        name = 'fr_dispersion_prop100';
    case 3
        objFun = @(eta,p) objFun_proportionalDissipation(FR({'eta',eta},p),0,false);
        name = 'fr_dispersion_prop0';
end
time = SSP_RK3;
export = struct('dat',false,'fig',true,'tikz',false);

%% Preprocess
try
    tblIn = readtable([name '.dat']);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(p);
tbl = array2table(zeros(I,7),'VariableNames',{'p','eta','c','exitFlag','relCost','A_T','CFL'});
parfor i = 1:I
    try
        options = optimset('PlotFcn',{
            @(varargin) plotFun_dispDiss(p(i),varargin{:})
        });
        [eta,badness,flag] = fminbnd(@(x) objFun(x,p(i)),-1,5,options);
        basis = FR({'eta',eta},p(i));
        badness = badness/objFun(0,p(i)) - 1; % relative to DG
        tbl{i,:} = [p(i) basis.eta basis.c flag badness basis.getOrder time.optimizeCFL(basis)];
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
        drawnow
        %%%
        if export.fig %#ok<PFBNS>
            saveas(gcf,[name '_' num2str(i) '.fig']);
        end
        if export.tikz
            cleanfigure
            matlab2tikz([name '_' num2str(i) '.tikz'],'showInfo',false)
        end
        %%%
    catch me
       warning("Run %d (p = %d) failed with error:\n '%s'",i,p(i),getReport(me))
    end
end

%% Postprocess
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[name '.dat'],'Delimiter','tab')
end