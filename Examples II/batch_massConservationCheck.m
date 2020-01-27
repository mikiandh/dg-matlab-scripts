clc
clear
%close all

% This script runs a batch of simulations to determine whether or not mass
% conservation (in the approximate solution, between initial and final
% times) is achieved.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Some settings
L = [-3,3];
% FUN = @(x) Functions.gauss(x,-3,3);
FUN = @(x) Functions.gauss(x,-1.5,0) + Functions.jump(x,0,1.5);
eqn = Advection;
time = @(limiter) SSP_RK3(0,12,eqn,limiter,.01,[]);
method = @(mesh) mesh.getSolutionMass;
% method = @(mesh) mesh.getSolutionNorm('mass');

%% Test matrix setup
%           {Method; Limiter; Nx}
data(:,1) = {DGIGA(10,1,0); []; 6};
data(:,2) = {DGIGA(10,2,0); []; 4};
data(:,3) = {DGIGA(10,2,1); []; 4};
data(:,4) = {DGIGA(10,3,0); []; 3};
data(:,5) = {DGIGA(10,3,1); []; 3};
data(:,6) = {DGIGA(10,3,2); []; 3};
data(:,7) = {DGIGA_AFC(10,1,0); Limiter.AFC(eqn); 6};
data(:,8) = {DGIGA_AFC(10,2,0); Limiter.AFC(eqn); 4};
data(:,9) = {DGIGA_AFC(10,2,1); Limiter.AFC(eqn); 4};
data(:,10) = {DGIGA_AFC(10,3,0); Limiter.AFC(eqn); 3};
data(:,11) = {DGIGA_AFC(10,3,1); Limiter.AFC(eqn); 3};
data(:,12) = {DGIGA_AFC(10,3,2); Limiter.AFC(eqn); 3};

%% Batch run
I = size(data,2);
data{end+3,end} = [];
parfor i = 1:I
% for i = 1:I
    data_i = data(:,i); % distribute
    mesh = Mesh(linspace(L(1),L(2),data_i{3}+1),data_i{1}.degree,data_i{1});
    data_i{1}.project(mesh,data_i{2},FUN);
    data_i{4} = method(mesh);
    solver = time(data_i{2});
    solver.launch(mesh,100,@(t,x) FUN(x));
    data_i{5} = method(mesh);
    data_i{6} = data_i{5} - data_i{4};
    data_i{1} = solver.plotData.space;
    data(:,i) = data_i; % gather
end

%% Table display
tab = cell2table(data');
tab.Properties.VariableNames = {'Method','Limiter','Nx','Start','End','Increase'};
disp(tab(:,[1 4:6]))