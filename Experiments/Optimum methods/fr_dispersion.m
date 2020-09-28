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
p = 5;%[2 3 4 5 7 10 14 19];
objFun = @(eta,p) objFun_Asthana2015(FR({'eta',eta},p));
time = SSP_RK3;
export = struct('name','fr_dispersion',...
    'dat',false,...
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
tbl = array2table(zeros(I,8),'VariableNames',{'p','eta','c','exitFlag','relCost','A_T','CFL','e1'});
parfor i = 1:I
    try
        options = optimset('PlotFcn',{
            @(varargin) plotFun_dispDiss(varargin{:},p(i),i,export)
        });
        [eta,badness,flag] = fminbnd(@(x) objFun(x,p(i)),-1,5,options);
        basis = FR({'eta',eta},p(i));
        badness = badness/objFun(0,p(i)) - 1; % relative to DG
        tbl{i,:} = [p(i) basis.eta basis.c flag badness basis.getOrder time.optimizeCFL(basis) basis.getCutoffWavenumber/basis.basisCount/pi];
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
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[export.name '.dat'],'Delimiter','tab')
end