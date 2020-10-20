clc
clear
close all

% This script generates a table comparing DGSEM with itself, at varying
% degrees, in terms of the common 'optimality criteria'.

%% Dependencies
addpath('../../Solver')
addpath('../../Basis')

%% Setup
p = [0 1 logspacei(2,19,8)];
time = SSP_RK3;
export = struct('name','dg_dispersion',...
    'dat',true,...
    'fig',true);

%% Preprocess
try
    tblIn = readtable([export.name '.dat']);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(p);
tbl = array2table(zeros(I,11),'VariableNames',{'p','Badness','relBadness','Order','relOrder','Resol','relResol','Cutoff','relCutoff','CFL_RK3','relCFL_RK3'});
parfor i = 1:I
    try
        optimum = DGSEM(p(i));
        tbl{i,:} = [
            p(i)
            objFun_Asthana2015(optimum)
            nan
            optimum.getOrder
            nan
            optimum.getResolvingWavenumber/optimum.basisCount/pi
            nan
            optimum.getCutoffWavenumber/optimum.basisCount/pi
            nan
            time.optimizeCFL(optimum)
            nan
        ]'; %#ok<PFBNS>
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
        %%%
        if export.fig %#ok<PFBNS>
            figure(i)
            optimum.displayDispDiss
            optimum.displayDispDissCombined
            saveas(gcf,[export.name '_' num2str(i) '.fig'])
            close hidden
        end
        %%%
    catch me
       warning("Run %d (p = %d) failed with error:\n '%s'",i,p(i),getReport(me))
    end
end

%% Postprocess
tbl{:,3:2:end} = tbl{:,2:2:end-1}./tbl{1,2:2:end-1} - 1; % relative changes over baseline
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[export.name '.dat'],'Delimiter','tab')
end