clc
clear
close all

% This script solves the following nonlinear, 3-dimensional, integer-valued
% optimization problem:
%
% Given DGIGA of 'J' basis functions, find the combination of 'k', 'p' and
% 's' parameters that minimizes some cost function associated with its
% modified dispersion and dissipation relations.

%% Dependencies
addpath('../../Solver')
addpath('../../Basis')

%% Setup
J = logspacei(2,19,8) + 1;
pMin = repelem(0,numel(J)); % minimum degree required
objFun = @objFun_Asthana2015;
name = 'dgiga_dispersion';
time = SSP_RK3;
export = struct('dat',true,'fig',false,'tikz',false);

%% Preprocess
try
    tblIn = readtable([name '.dat']);
    J(any(J == tblIn.J & pMin == tblIn.pMin)) = []; % drop repeated runs
    J(J < 2) = []; % drop FV methods
catch
    tblIn = table; % empty table
end

%% Minimization
I = numel(J);
tbl = array2table(zeros(I,26),'VariableNames',{
    'J'
    'pMin'
    'k'
    'p'
    's'
    'Cond'
    'relCondBern'
    'relCondLagr'
    'Goodness'
    'relGoodnessBern'
    'relGoodnessLagr'
    'DispDiss'
    'relDispDissBern'
    'relDispDissLagr'
    'Order'
    'relOrderBern'
    'relOrderLagr'
    'Resol'
    'relResolBern'
    'relResolLagr'
    'Cutoff'
    'relCutoffBern'
    'relCutoffLagr'
    'CFL_RK3'
    'relCFLRK3Bern'
    'relCFLRK3Lagr'
    });
% In addition, store ALL candidate scores
%%%tblExtra = array2table(zeros(I,4),'VariableNames',{'f','k','p','s'});
F = cell(I,1); % array of candidates
for i = 1:I
    try
        % Generate all possible candidates:
        [p,s] = meshgrid(1:J(i)-1,0:J(i)-2);
        k = (J(i)-s-1)./(p-s);
        isOk = s < p & floor(k) == k & p >= pMin(i) & (k > 1 | s == p-1);
        k(~isOk) = [];
        p(~isOk) = [];
        s(~isOk) = [];
        % Evaluate cost function for each:
        f = zeros(1,sum(isOk(:)));
        %%%
        figure
        xlabel('\kappa')
        yyaxis left
        plot([0 pi*J(i)],[0 pi*J(i)],'--','DisplayName','Exact (disp.)')
        ylabel('\Re(\omega)')
        yyaxis right
        plot([0 pi*J(i)],[0 0],'--','DisplayName','Exact (diss.)')
        ylabel('\Im(\omega)')
        hold on
        %%%
        % Compute (in parallel):
        D = parallel.pool.DataQueue;
        D.afterEach(@plotFun_single);
        condEstNums = zeros(numel(f),1);
        parfor j = 1:numel(f)
            basis = DGIGA(k(j),p(j),s(j));
            f(j) = objFun(basis);
            condEstNums(j) = condest(basis.massMatrix);
            %%%
            [w0,w] = basis.getCombinedWavenumbers;
            send(D, {j,w0(w0 >= 0),w(w0 >= 0),basis.getName})
            %%%
        end
        % Store them all:
        F{i} = table(repmat(J(i),numel(f),1),k.',p.',s.',condEstNums,f.','VariableNames',{'J','k','p','s','Condest','Badness'});
        % Extract optimum:
        [~,j] = min(f);
        %%%
        yyaxis left
        set(findobj(get(gca,'Children'),'Tag',['disp_' num2str(j)]),'LineStyle','-');
        yyaxis right
        set(findobj(get(gca,'Children'),'Tag',['diss_' num2str(j)]),'LineStyle','-');
        title(sprintf('DGIGA, J = %d (%d combinations)',J(i),numel(f)))
        %%%
        optimum = DGIGA(k(j),p(j),s(j));
        baselines = [DGIGA(1,p(k == 1)) DGSEM(p(k == 1))];
        [~,~,dispDiss0] = optimum.getDispDissRatio;
        [~,~,dispDiss1] = baselines(1).getDispDissRatio;
        [~,~,dispDiss2] = baselines(2).getDispDissRatio;
        tbl{i,:} = [
            J(i)
            pMin(i)
            k(j)
            p(j)
            s(j)
            condEstNums(j)
            condEstNums(k == 1)
            condest(baselines(2).massMatrix)
            -f(j)
            -f(k == 1)
            -objFun_Asthana2015(baselines(2))
            dispDiss0
            dispDiss1
            dispDiss2
            optimum.getOrder
            baselines(1).getOrder
            baselines(2).getOrder
            optimum.getResolvingWavenumber/optimum.basisCount/pi
            baselines(1).getResolvingWavenumber/baselines(1).basisCount/pi
            baselines(2).getResolvingWavenumber/baselines(2).basisCount/pi
            optimum.getCutoffWavenumber/optimum.basisCount/pi
            baselines(1).getCutoffWavenumber/baselines(1).basisCount/pi
            baselines(2).getCutoffWavenumber/baselines(2).basisCount/pi
            time.optimizeCFL(optimum)
            time.optimizeCFL(baselines(1))
            time.optimizeCFL(baselines(2))
        ]';
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
        drawnow
        %%%
        if export.fig
            saveas(gcf,[name '_' num2str(i) '.fig']);
        end
        if export.tikz
            cleanfigure
            matlab2tikz([name '_' num2str(i) '.tikz'],'showInfo',false)
        end
        %%%
    catch me
       warning("Run %d failed with error:\n '%s'",i,getReport(me))
    end
end

%% Postprocess
if ~isempty(tbl)
    tbl{:,[7:3:end 8:3:end]} = (tbl{:,[6:3:end-2 6:3:end-2]} - tbl{:,[7:3:end 8:3:end]})./abs(tbl{:,[7:3:end 8:3:end]}); % relative changes over baseline
end
tblOut = sortrows([tblIn; tbl],'J');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[name '.dat'],'Delimiter','tab')
    writetable(sortrows(vertcat(F{:}),'Badness'),[name '_all' '.dat'],'Delimiter','tab');
end