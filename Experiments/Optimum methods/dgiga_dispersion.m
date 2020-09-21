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
J = [3 4 5 6]; % 8 11 15 20];
pMin = 0*J; % minimum degree required
switch 2
    case 1
        objFun = @(basis) objFun_peakPosition(basis,.5);
        name = 'dgiga_dispersion_peak50';
    case 2
        objFun = @(basis) objFun_proportionalDissipation(basis,1,true);
        name = 'dgiga_dispersion_prop100';
end
time = SSP_RK3;
export = struct('dat',true,'fig',true,'tikz',false);

%% Preprocess
try
    tblIn = readtable([name '.dat']);
    J(any(J == tblIn.J & pMin == tblIn.pMin)) = []; % drop repeated runs
    J(J < 2) = []; % drop FV methods
catch
    tblIn = table; % empty table
end
tbl = array2table(zeros(numel(J),8),'VariableNames',{'J','pMin','k','p','s','relCost','A_T','CFL'});

%% Minimization
I = numel(J);
parfor i = 1:I
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
        hold on
        xlabel('\kappa')
        %%%
        for j = 1:numel(f)
            basis = DGIGA(k(j),p(j),s(j));
            f(j) = objFun(basis);
            %%%
            [w,w0] = basis.getFourierFootprint;
            w = 1i*w(1,w0 > 0);
            w0(w0 < 0) = [];
            yyaxis left
            plot(w0,real(w),':','DisplayName',sprintf('f = %g',f(j)),'Tag',['disp_' num2str(j)]);
            ylabel('\Re(\omega)')
            yyaxis right
            plot(w0,imag(w),':','Tag',['diss_' num2str(j)]);
            ylabel('\Im(\omega)')
            %%%
        end
        % Extract optimum:
        [~,j] = min(f);
        basis = DGIGA(k(j),p(j),s(j));
        %%%
        yyaxis left
        set(findobj(get(gca,'Children'),'Tag',['disp_' num2str(j)]),'LineStyle','-');
        set(refline(1,0),'LineStyle','--')
        yyaxis right
        set(findobj(get(gca,'Children'),'Tag',['diss_' num2str(j)]),'LineStyle','-');
        set(refline(0,0),'LineStyle','--')
        title(sprintf('DGIGA, J = %d (%d combinations)',J(i),numel(f)))
        %%%
        tbl{i,:} = [J(i) pMin(i) k(j) p(j) s(j) f(j)./f(k==1)-1 basis.getOrder time.optimizeCFL(basis)]; %#ok<PFBNS>
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
       warning("Run %d failed with error:\n '%s'",i,getReport(me))
    end
end

%% Postprocess
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[name '.dat'],'Delimiter','tab')
end