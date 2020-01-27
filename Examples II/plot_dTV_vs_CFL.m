clc
clear
close all

% This script generates a plot of CFL number vs. change in Total Variation
% for a fixed set of remaining parameters.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Math')
addpath('../Grid')

%% Parameters
N = 20;
CFLs = linspace(.5,5,N);
Nk = 10; % number of patches
p = 1; % degree of approximation
method = DGIGA(2);
timeIntegrator = SSP_RK4_10(0,[],50,Advection,[]);

%% Setup
fun = @(t,x) Functions.gauss(x);
mesh = Mesh(linspace(0,1,Nk+1),p,method);

%% Data sampling (parallel)
tic
parfor i = 1:N
    dTVs(:,i) = dTV_function(method,timeIntegrator,mesh,CFLs(i),fun);
    fprintf(1,'Run %d out of %d (worker %d)\n',i,N,get(getCurrentTask(), 'ID'))
end
toc

%% dTV vs. CFL
figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(CFLs,dTVs,'-+b');
xlabel('CFL')
ylabel('\DeltaTV')
% id = find(ischange(dTVs,'variance','MaxNumChanges',1));
% zoomPlot(CFLs,dTVs,[CFLs(max(1,id-50)) CFLs(id)],[.2 .33 .33 .33],[1 2 3 4]);
ylim([-inf 5])

%% abs(dTV) vs. CFL
figure('Renderer', 'painters', 'Position', [1100 600 800 400])
semilogy(CFLs,abs(dTVs),'-+r');
xlabel('CFL')
ylabel('| \DeltaTV |')

%% Auxiliary functions
function dtv = dTV_function(method,timeIntegrator,mesh,CFL,fun)
    % Total variation difference
    if method.isHybrid && mesh.maxDegree > 0
        Ns = method.nonzeroSpanCount;
    else
        Ns = 1;
    end
    dt = CFL/mesh.elementCount/Ns*2.89/(mesh.maxDegree^2 + 3*mesh.maxDegree + 2.89);
    method.project(mesh,[],@(x) fun(0,x));
    dtv = mesh.getTotalVariation;
    timeIntegrator.timeNow = 0;
    if timeIntegrator.launchFixedTimeStep(dt,mesh,500,fun)
        dtv = mesh.getTotalVariation - dtv;
    else
        dtv = 1e300;
    end
    close
end