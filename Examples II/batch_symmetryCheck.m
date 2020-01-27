clc
clear
%close all

% This script runs a batch of advection simulations to determine wether or
% not each discretization fullfils any or both of the following symmetry 
% measures around the origin:
%  1) Mirror symmetry (even function)
%  2) Analogue symmetry (upwind and downwind sides match, respectively, 
%     after switching the advection direction)
%
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
FUN = @(x) Functions.gauss(x,-3,3);

%% Test matrix setup
%           {Method; Limiter; Nx}
data(:,1) = {DG(0); 'none'; 48};
data(:,2) = {DGIGA(5,3,0); 'none'; 3};
data(:,3) = {DGIGA(5,3,1); 'none'; 4};
data(:,4) = {DGIGA(5,3,2); 'none'; 6};

%% Batch run
I = size(data,2);
data{end+2,end} = []; % type 1 symmetry error, type 2 symmetry error
% parfor i = 1:I
for i = 1:I
    data_i = data(:,i); % distribute
    mesh = Mesh(linspace(L(1),L(2),data_i{3}+1),data_i{1}.degree,data_i{1});
    for a = [1 -1]
        eqn = Advection(a);
        limiter = Limiter.(data_i{2})(eqn);
        data_i{1}.project(mesh,limiter,FUN);
        time = SSP_RK3(0,12,eqn,limiter,.05,[]);
        time.launch(mesh,0,@(t,x) FUN(x));
        switch a
            case 1
                data_i{4} = mesh.getErrorNorm(@(x) mesh.sample(-x));
                aux = mesh;
            otherwise
                data_i{5} = mesh.getErrorNorm(@(x) aux.sample(x));
        end
    end
    data_i{1} = time.plotData.space;
    data(:,i) = data_i; % gather
end

%% Table display
tab = cell2table(data');
tab.Properties.VariableNames = {'Method','Limiter','Nx','Symmetry','Analogy'};
disp(tab(:,[1 4:5]))