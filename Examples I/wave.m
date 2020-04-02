clc
clear
%close all
%path(pathdef)

% This script solves the wave equation in 1D.

%% Dependencies
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Basis')
addpath('../Mesh')
addpath('../Math')
addpath('../Extra')

%% Discretization (equation + solution + domain)
mesh = Mesh(DG(3),[0 1],[Reflective(@(t) .5*sin(2*pi*t)) Reflective],50);

%% Solver
norms = Norm(["Mass","L2"]);
solver = SSP_RK3(Wave,[0 4],...
    'courantNumber',.1,...
    'limiter',BSB('Sensor',KXRCF),...
    'showSensor',true,'showLimiter',true,'norms',norms,...
    'exactSolution',@constantIC,'iterSkip',50);

%% Time-integration
solver.initialize(mesh)
solver.launch(mesh);

%% Postprocessing
solver.physics.displayData(num2cell(norms),{'Norm','tEnd'},[norms.vals])

%% Auxiliary functions
function y = constantIC(~,x)
y = zeros(2,length(x));
%%% y(2,:) = 1;
end

function y = pickIC(~,x)
% Initial condition that models a smooth displacement pulse.
y = zeros(2,length(x));
y(1,:) = exp(-32*(x).^2);
end

function y = hammerIC(~,x)
% Initial condition that models a sudden perturbation in displacement rate 
% on a portion of an initially still medium.
y = zeros(2,length(x));
y(2,:) = heaviside(x-2/3+1) - heaviside(x-2*2/3+1);
end

function y = standingIC(~,x)
% Initial condition that will generate a smooth standing wave.
y = zeros(2,length(x));
y(2,:) = sin(2*x*pi);
end

function y = combinedIC(~,x)
% Initial condition that models the combination of a 'pick' (smooth pulse
% in displacement) and a 'hammer strike' (sudden change in displacement 
% rate) over non-overlapping regions of the domain.
y = zeros(2,length(x));
y(1,:) = Functions.gauss(x,0,.5);
y(2,:) = Functions.jump(x,.5,1);
end