clc
clear
%close all
%path(pathdef)

% This script computes the L2 error of the DGSEM solution to the Advection
% equation with a monochromatic signal of various wavenumbers as inital
% condition. The objective is to experimentally verify the spectral
% properties measured via MWA by comparing L2 errors.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Settings
p = [0 1 3 5];
Nx = 16;
dt = 1e-3;
tEnd = .1;

%% Test tuple setup
data = assembleTestTuple(...
    {DGSEM}... space (2)
    ,{1 3 5}... degree (3)
    ,{16}... #patches (4)
    ,{@SSP_RK3}... time (5)
    ,{1e-4}... time-step size (6)
    ,{.05}... final time (7)
    ,{nan}... wavenumbers (8)
    ,{nan}... L2 errors (9)
    ,{nan}... WCT (10)
    );

%% Parallel batch run
I = size(data,1);
parfor i = 1:I
% for i = 1:I
    data_i = data(i,:);
    mesh = Mesh(linspace(0,2*pi,data_i{4}+1),data_i{3},data_i{2});
    time = data_i{5}(0,[],data_i{7},Advection,[]);
    tic
    [data_i{8}, data_i{9}] = computeErrors(mesh,time,data_i{6});
    data_i{10} = toc;
    % Gather data back:
    data(i,:) = data_i;
end

%% Combined plot
c = distinguishable_colors(I,{'w','k'});
figure('Renderer', 'painters', 'Position', [1000 100 700 400])
legend('-DynamicLegend','Location','Best')
hold all
for i = 1:I
    name = sprintf('%s; p = %d, N_x = %d',...
        class(data{i,2}),data{i,3},data{i,4});
    plot(data{i,8},data{i,9},'*-','Color',c(i,:),'DisplayName',name);
end
hold off
ylabel('L2 error')
xlabel('$$\kappa \frac{\Delta x}{p + N_\sigma}$$','Interpreter','latex')

%% Auxiliary functions
function [kHat,errs] = computeErrors(mesh,time,dt)
k = 1:1:mesh.dofCount/2;
if mod(k(end),2)
    error('#DOFs must be even. Try different settings.')
end
kHat = k*2*pi./mesh.dofCount;
errs = nan(size(kHat));
for i = 1:length(k)
    time.timeNow = 0;
    mesh.bases.project(mesh,[],@(x) sin(k(i)*x));
    time.launchFixedTimeStep(dt,mesh,0,@(t,x) sin(k(i)*(x-t)));
    errs(i) = mesh.getErrorNorm(@(x) sin(k(i)*(x-time.timeStop)));
end
end