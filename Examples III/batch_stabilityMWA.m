clc
clear
%close all
%path(pathdef)

% This script approximates the maximum time-step size bounds for given time
% and space semi-discretization combinations. Procedure adapted from van
% den Abeele 2009 (PhD thesis), section B3.
%
% The stability limit occurs when one point in the scaled Fourier footprint
% of the spatial scheme reaches the edge of the time scheme's stability 
% region - i.e. the amplification factor of the full (space and time)
% discretization first reaches unity as the time-step size decreases from 
% infinity.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Solver')

%% Test tuple setup
data = assembleTestTuple(...
    {SSP_RK1 SSP_RK2 SSP_RK3}... time semi-discretization scheme (2)
    ,{DGIGA(1)}... space semi-discretization scheme (3)
    ,{1}... upwinding ratio (4)
    ,{0 1 2}... degree (5)
    ,{60}... #patches (6)
    ,{nan}... #spans (7)
    ,{nan}... #DOFs (8)
    ,{10}... critical Courant number (9)
    ,{nan}... root-finder exit flag (10)
    ,{nan}... WCT (11)
    );

%% Batch run
I = size(data,1);
fprintf(1,...
    '%7s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%12s\n',...
    '#/Total','Time','Space','Upwind','p','#patches','#spans','#DOFs',...
    'CFL','Success?','WCT(s)');
% parfor i = 1:I
for i = 1:I
    data_i = data(i,:);
    if data_i{3}.isHybrid
        data_i{7} = data_i{3}.nonzeroSpanCount;
    else
        data_i{7} = 1;
    end
    data_i{8} = data_i{6}*(data_i{5}+data_i{7});
    theta = -1i*MWA_eigen_full(data_i{[6 5 3 4]}); % spatial eigenvalues
    tic
    [data_i{[9 10]}] =...
        optimizeCFL(data_i{9},theta(:),data_i{2}.amplificationFactorFun);
    data_i{11} = toc;
    % Print solution:
    fprintf(1,...
    '%3d/%-3d\t%8s\t%8s\t%8f\t%8d\t%8d\t%8d\t%8d\t%8.4g\t%8d\t%12g\n',...
    data_i{1},I,class(data_i{2}),class(data_i{3}),data_i{4:end});
    % Gather data back:
    data(i,:) = data_i;
end

%% Table display
clc
disp(cell2table(data,'VariableNames',...
    {'i','time','space','upwind','p','Nx','Ns','N',...
    'CFL','exitFlag','WCT'}))

%% Nonlinear constrained optimization subroutine
function [CFL,exitFlag] = optimizeCFL(CFL,theta,G)
% This function sets up a constrained optimization problem, and solves it 
% via fmincon. The goal is to determine the CFL coefficient that maximizes
% the L1 norm of the amplification factors - all eigenmodes and wavenumbers
% considered.
%
% Nonlinear function to minimize:
    function norm = fun(CFL)
        norm = - sum(abs(G(theta*CFL)));
    end
% Nonlinear constraint:
    function [c,ceq] = nonlcon(CFL)
        c = max(abs(G(theta*CFL))) - 1;
        ceq = [];
    end
% Solve the problem:
[CFL,~,exitFlag] = fmincon(@fun,CFL,[],[],[],[],0,[],@nonlcon,...
    optimoptions('fmincon','Display','off'));
end