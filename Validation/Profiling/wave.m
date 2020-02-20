clc
clear
close all
%path(pathdef)

% This script compares wall clock timings for various discretizations when 
% applied to the wave equation test case.

%% Dependencies
addpath('../../Extra')
addpath('../../Discretization')
addpath('../../Limiters')
addpath('../../Physics')
addpath('../../Solver')
addpath('../../Grid')
addpath('../../Math')

%% Parameters
params = {...
    DG(1), 10
    DG(1), 20
    DG(2), 10
    DG(2), 20};

%% Common set up
CFL = .1; % Courant number
eqn = Wave;
fun = @(x) combinedIC(x);
limiter = [];
WCTs = nan(size(params,1),1);

%% Solution
for i = 1:length(WCTs)
    Ne = params{i,2}+1;
    method = params{i,1};
    p = method.degree;
    mesh = Mesh(linspace(0,1,Ne+1),p,method);
    method.project(mesh,limiter,fun);
    timeIntegrator = SSP_RK3(0,2,eqn,limiter,CFL,[]);
    tic
    timeIntegrator.launch(mesh,0,@(t,x) fun(x));
    WCTs(i) = toc;
end

%% Results
results = [cell2table(params) table(WCTs)];
disp(results)

%% Functions
function y = combinedIC(x)
% Initial condition that models the combination of a 'pick' (smooth pulse
% in displacement) and a 'hammer strike' (sudden change in displacement 
% rate) over non-overlapping regions of the domain.
y = zeros(2,length(x));
y(1,:) = Functions.gauss(x,0,.5);
y(2,:) = Functions.jump(x,.5,1);
end