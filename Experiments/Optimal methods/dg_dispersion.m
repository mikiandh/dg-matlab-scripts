clc
clear

% This script generates a table comparing DGSEM with itself, at varying
% degrees, in terms of the common 'optimality criteria'.

%% Setup
p = [0 1 logspacei(2,19,8)];
time = SSP_RK3;
export = struct('name','dg_dispersion',...
    'dat',true,...
    'fig',false);

%% Preprocess
try
    tblIn = readtable([export.name '.dat']);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(p);
tbl = array2table(zeros(I,13),'VariableNames',{
    'p'
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
        optimum = DGSEM(p(i));
        [~,~,R] = optimum.getDispDissRatio;
        tbl{i,:} = [
            p(i)
            -objFun_Asthana2015(optimum)
            nan
            R
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This plots and saves (if active) the dispersion and dissipation
        % relations obtained from using simple MWA (physical mode), and
        % combined-mode analysis. Neither are actually representative
        % of actual behaviour; only phase shift and amplification factor
        % are.
        if export.fig %#ok<PFBNS>
            figure(i)
            optimum.displayDispDiss
            optimum.displayDispDissCombined
            saveas(gcf,[export.name '_' num2str(i) '.fig'])
            close hidden
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch me
       warning("Run %d (p = %d) failed with error:\n '%s'",i,p(i),getReport(me))
    end
end

%% Postprocess
if ~isempty(tbl)
    tbl{:,3:2:end} = (tbl{:,2:2:end-1} - tbl{1,2:2:end-1})./abs(tbl{1,2:2:end-1}); % relative changes over baseline
end
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[export.name '.dat'],'Delimiter','tab')
end