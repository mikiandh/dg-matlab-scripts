clc
clear
%close all

% This script solves the (inviscid) Burgers' equation.

%% Dependencies
addpath('../Limiting')
addpath('../Physics')
addpath('../Solver')
addpath('../Basis')
addpath('../Mesh')
addpath('../Math')
addpath('../Extra')

%% Discretization
mesh = Mesh(DGIGA(5,2),[-1 1],[Farfield(2) Transmissive],20);

%% Solver
solver = SSP_RK3(Burgers,[0 .4],...
    'courantNumber',.01,...
    'exactSolution',@hesthaven,'iterSkip',60);

%% Initial condition
solver.initialize(mesh)

%% Time-integration
solver.launch(mesh);

%% Initial conditions
function y = riemann1A(~,x) % right-going shock
y = 1 - 0.8*heaviside(x);
end

function y = riemann1B(~,x) % right-going expansion
y = 0.2 + 0.8*heaviside(x);
end

function y = riemann2A(~,x) % left-going expansion
y = -1 + 0.8*heaviside(x);
end

function y = riemann2B(~,x) % left-going shock
y = 1 - 1.2*heaviside(x);
end

function y = riemann3A(~,x) % shock-shock collision (right-going wins)
y = 1 - 1.2*heaviside(x);
end

function y = riemann3B(~,x) % shock-shock collision (left-going wins)
y = 0.2 - heaviside(x);
end

function y = riemann3C(~,x) % shock-shock collision (standing shock)
y = 1 - 2*heaviside(x);
end

function y = riemann4A(~,x) % transonic expansion (right-going)
y = -0.2 + heaviside(x);
end

function y = riemann4B(~,x) % transonic expansion (left-going)
y = -0.8 + heaviside(x);
end

function y = gaussIC(~,x)
y = Functions.gauss(x,-1,1);
end

function y = jumpIC(~,x)
y = Functions.jump(x,-1,1);
end

function y = halfJumpIC(~,x)
y = x.*(heaviside(x) - heaviside(x-2));
end

function y = sineIC(~,x)
y = sin(x) + 0.5;
end

function y = combinedIC(~,x) % perfect initial condition; L = [-1,2], tEnd = 2
y = Functions.gauss(x,-1,0) + Functions.jump(x,0,1);
end

function y = toroIC(~,x) % from Toro, 2009 (p. 196); L = [0 1.5], tEnd = .5
y = -.5 + 1.5*heaviside(x-.5) - 1*heaviside(x-1);
end

function y = leveque1(t,x) % Leveque (expansion wave), 2002 (p. 231); L = [-3 3], tEnd = .5.
% NOTE: cancellation of fluxes may occur if -qL = qR
% NOTE: IGA breaks whenever N is even if shock is at x == 0
if t == 0
    y = -1 + 3*heaviside(x);
else
    y = -1 + (1+x/t).*heaviside(x+t) + (2-x/t).*heaviside(x-2*t);
end
end

function y = leveque2(t,x) % Adapted from Leveque, 2002 (p. 238); L = [-1 15], tEnd = 14.
y = 2 - 2*heaviside(x-t);
end

function y = leveque3(~,x) % Inspired by Leveque, 2002 (p. 223); L = [-8 8], tEnd = 6.
y = Functions.gauss(x,-4,4).*Functions.tones(x,[1 4 8],[1 2 3],-8,8);
end

function y = hesthaven(t,x) % Hesthaven & Warburton, 2008 (p. 141); L = [-1 1], tEnd = 0.8.
y = 2 - heaviside(x+.5-1.5*t);
end