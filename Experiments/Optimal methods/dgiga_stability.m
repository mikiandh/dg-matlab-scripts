clc
clear

% This script solves the following nonlinear, 3-dimensional, integer-valued
% optimization problem:
%
% Given DGIGA of 'J' basis functions, find the combination of 'k', 'p' and
% 's' parameters that maximizes the maximum allowable CFL number.

%% Setup
J = [3 4 5 6]; % 8 11 15 20];
pMin = 0*J; % minimum degree required
time = SSP_RK3;
name = 'dgiga_stability';
export = struct('dat',true,'fig',false,'tikz',false);

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
        % Evaluate fitness function for each:
        CFL = zeros(1,sum(isOk(:)));
        for j = 1:numel(CFL)
            CFL(j) = time.optimizeCFL(DGIGA(k(j),p(j),s(j))); %#ok<PFBNS>
        end
        % Extract optimum:
        [~,j] = max(CFL);
        tbl{i,:} = [J(i) pMin(i) k(j) p(j) s(j) -CFL(j)./CFL(k==1)+1 DGIGA(k(j),p(j),s(j)).getOrder CFL(j)];
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
    catch me
       warning("Run %d failed with error:\n '%s'",i,getReport(me))
    end
end

%% Postprocess
tblOut = sortrows([tblIn; tbl],'J');
clc
disp(tblOut)
if export.dat
    writetable(tblOut,[name '.dat'],'Delimiter','tab')
end