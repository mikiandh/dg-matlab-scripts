clc
clear
close all

% This script runs a batch of simulations and saves their parameters and 
% corresponding L2 error as a table in a file.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')

%% Batch parameters
filename = 'gauss_DGIGA_RK3.dat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 dimensions, no more no less, and in this order %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = num2cell(logspace(-4,0,40)); %num2cell([1e-3 reshape((2:10)'.*10.^(-3:-2),1,[])]); % time-step size
Nk = {1}; % number of patches
Nl = {10}; % number of knot spans
p = num2cell(1:29); % degree of approximation

%% Batch run setup
data = assembleTestMatrix(dt,Nk,Nl,p,{nan},{nan},{nan}); % row: simulation; column: variable

%% Batch simulation runs
eqn = Advection(1);
FUN = @(t,x) exp(-18*(.5*sqrt(2*pi))^2*(x-.5).^2);
tEnd = 1;
fprintf(1,'%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n','#','dt','#patches','#spans','p','#DOFs','Error','WCT(s)');
for runId = 1:size(data,1)
    tic
    % /////////////////////////////////////////////////////////////////////
    method = DGIGA(data(runId,3)); % discretization
    limiter = []; % Limiter.AFC(eqn);
    % ---------------------------------------------------------------------
    %method = DGSEM;
    %limiter = []; %Limiter.Krivodonova(eqn);
    % /////////////////////////////////////////////////////////////////////
    mesh = Mesh(linspace(0,1,data(runId,2)+1),data(runId,4),method); % grid
    method.project(mesh,limiter,@(x) FUN(0,x)); % initial solution
    timeIntegrator = SSP_RK3(0,[],tEnd,eqn,limiter); % solver
    timeIntegrator.launchFixedTimeStep(data(runId,1),mesh,1000,FUN); % launch run
    data(runId,5) = mesh.dofCount;
    data(runId,6) = mesh.getErrorNorm(@(x) FUN(tEnd,x)); % L2 error norm
    data(runId,7) = toc;
    fprintf(1,'%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n',runId,data(runId,:));
    close
end

%% Save data to disk
if ~isempty(filename)
    fileID = fopen(filename,'w');
    fprintf(fileID,'%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n','#','dt','#patches','#spans','p','#DOFs','Error','WCT(s)');
    for runId = 1:size(data,1)
        fprintf(fileID,'%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n',runId,data(runId,:));
    end
    fclose(fileID);
end

%% Plot data
figure()
subplot(2,1,1)
loglog(data(1:end,1),data(1:end,6),'*-b')
subplot(2,1,2)
loglog(data(1:end,1),abs(gradient(data(1:end,6))),'*-r')

%% Auxiliary functions
function data = assembleTestMatrix(varargin)
    % Assembles all combinations of an arbitrary number of cell arrays of 
    % any number of variables each into a test matrix, such that each row 
    % is a combination and each column is an input cell array.
    %
    % Based on the "fundamental principle of counting": 
    % https://stackoverflow.com/questions/3536833/arbitrary-number-of-nested-loops
    %
    % Preallocation:
    N = cellfun('length',varargin); % number of data points in each category
    M = prod(N); % total number of combinations taking one data point from each category
    data = zeros(M,nargin);
    % Loop over combinations:
    n = ones(size(N)); % set each category's counter
    for i = 1:M
        % Loop over categories:
        for j = 1:nargin
            data(i,j) = varargin{j}{n(j)}; % select one cell from each category
        end
        % Advance each category's counter:
        for k = 1:nargin
            n(k) = n(k) + 1;
            if n(k) > N(k)
                n(k) = 1; % if it surpasses its maximum, reset it and advance next category
            else
                break % if not, do not advance the remaining categories
            end
        end
    end
end