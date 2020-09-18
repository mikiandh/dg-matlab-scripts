%clc
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
p = 5;
time = SSP_RK3;
filename = 'fr_dispersion.dat';

%% Preprocess
try
    tblIn = readtable(filename);
    p(any(p == tblIn.p)) = []; % drop repeated degrees
catch
    tblIn = table; % empty table
end
tbl = array2table(zeros(numel(p),5),'VariableNames',{'p','eta','c','badness','exitFlag'});

%% Minimization
I = numel(p);
for i = 1:I
    try
        options = optimset('PlotFcn',{
            @(varargin) plotDispDiss(p(i),varargin{:})
%             @(varargin) plotEigenvalues(p(i),varargin{:})
        });
        [eta,badness,flag] = fminbnd(@(x) objFun2(x,p(i)),0,5,options);
        basis = FR({'eta',eta},p(i));
        tbl{i,:} = [p(i) basis.eta basis.c badness flag];
        fprintf('\nRun %d of %d:\n',i,I)
        disp(tbl(i,:))
    catch
       warning('Run %d (p = %d) failed.',i,p(i))
    end
end

%% Postprocess
return
tblOut = sortrows([tblIn; tbl],'p');
clc
disp(tblOut)
writetable(tblOut,filename,'Delimiter','tab')

%% Objective functions
function badness = objFun1(eta,p)
% Minimize over/undershoot and push to the right the dispersion peak.
%
r = .5; % overshoot/range weighting (0 to 1)
%
[z,k] = FR({'eta',eta},p).getFourierFootprint;
w = 1i*z(1,k > 0);
k(k <= 0) = [];
[y,id] = max(real(w));
badness = r*abs(y - k(id)) + (1-r)*(k(end) - k(id));
end

function f = objFun2(eta,p)
% Minimize dispersion error but encourage dissipation in proportion to the
% former.
%
r = 1; % ratio between dissipation and dispersion errors (after cutoff)
%
[z,k] = FR({'eta',eta},p).getFourierFootprint;
w = 1i*z(1,k > 0);
k(k <= 0) = [];
err = real(w) - k;
nc = find(err >= 0,1,'last'); % cutoff wavemode (right-most intercept)
w0 = k - r*1i*(k > k(nc)).*abs(err); % target Block wave
f = norm(w - w0); % badness
end

%% Plot functions
function stop = plotDispDiss(p,eta,optimValues,state,varargin)
% Plots the dispersion and dissipation relations of the current iteration.
% Compares them with those of the previous iteration, DG's, and exact ones.
stop = false;
switch state
    case 'init' % set up the plot
        [z,k] = FR('DG',p).getFourierFootprint;
        hold on
        yyaxis left
        plot(k,k,'--k','DisplayName','Exact','Tag','exactDisp')
        plot(k,-imag(z(1,:)),':k','DisplayName','eta = 0','Tag','refDisp')
        ylabel \Re(\omega*)
        yyaxis right
        plot(k,0*k,'--k','Tag','exactDiss')
        plot(k,real(z(1,:)),':k','Tag','refDiss')
        ylabel \Im(\omega*)
        xlabel \kappa*
        xlim([0 inf])
        title(sprintf('Degree: %d',p))
    case 'iter'
        [z,k] = FR({'eta',eta},p).getFourierFootprint;
        if optimValues.iteration > 0
            for side = ["left" "right"]
                yyaxis(side)
                delete(findobj(get(gca,'Children'),'Tag','old'));
                hNew = findobj(get(gca,'Children'),'Tag','new');
                hOld = copy(hNew,gca);
                set(hOld,{'Color','Tag'},{'b','old'});
                delete(hNew)
            end
        end
        yyaxis left
        plot(k,-imag(z(1,:)),'-r',...
            'DisplayName',sprintf('eta = %g (f = %g)',eta,optimValues.fval),...
            'Tag','new')
        yyaxis right
        plot(k,real(z(1,:)),'-r','Tag','new')
        yyaxis left
        legend(get(gca,'Children'),'Location','South')
end
end

function stop = plotEigenvalues(p,eta,optimValues,state,varargin)
% Plots the fourier footprint of the current iteration.
stop = false;
switch state
    case 'init' % set up the plot
        z = FR('DG',p).getFourierFootprint;
        set(gca,{'XAxisLocation','YAxisLocation'},{'origin','origin'})
        plot(z(:),'.k','DisplayName','eta = 0');
        hold on
        title(sprintf('Degree: %d',p))
        xlabel \Re(z*)
        ylabel \Im(z*)
        legend('-DynamicLegend','Location','SouthWest')
    case 'iter'
        z = FR({'eta',eta},p).getFourierFootprint;
        if optimValues.iteration > 0
            delete(findobj(get(gca,'Children'),'Tag','old'));
            hNew = findobj(get(gca,'Children'),'Tag','new');
            hOld = copy(hNew,gca);
            set(hOld,{'Color','Tag'},{'b','old'});
            delete(hNew)
        end
        plot(z(:),'.r','DisplayName',sprintf('eta = %g (f = %g)',eta,optimValues.fval),'Tag','new');
end
end