clc
clear
close all

% This script runs a batch of simulations and determines the maximum 
% allowable CFL number for every combination of input parameters. Results 
% are written to a file.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Math')
addpath('../Grid')

%% Batch run setup
filename = 'CFL.dat';
data = assembleTestTuple(...
    {SSP_RK4_10(0,[],50,Advection,[])}... temporal scheme (2)
    ,{DGIGA(2)}... spatial scheme (3)
    ,{1}... polynomial degree (4)
    ,{nan}... #patches, preallocation (5)
    ,{nan}... #spans, preallocation (6)
    ,{30}... degrees of freedom, target (7)
    ,{'fzero'}... root/minimum search algorithm (8)
    ,{1e-16}... CFL tolerance (9)
    ,{[4 9]}... CFL search interval (10)
    ,{nan}... dt, preallocation (11)
    ,{nan}... dTV, preallocation (12)
    ,{nan}... exit status flag, preallocation (13)
    ,{nan}... WCT, preallocation (14)
    ,{nan}... CFL, preallocation (15)
    );

%% Parallel batch run
I = size(data,1);
fprintf('Launching %d simulations, using %d workers.\n',I,maxNumCompThreads);
varNames = {'#','Total','Time','Space','Degree','Patches','Spans','DOFs','Solver','Tol','Range','dt','dTV','Status','WCT(s)','CFL'};
varNamesSpecs = '%3s/%-4s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-14s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\n';
varSpecs = '%3d/%-4d\t%-8s\t%-8s\t%-8d\t%-8d\t%-8d\t%-8d\t%-8s\t%-8.6g\t[%5.4g, %-5.4g]\t%-8.6g\t%-8.4g\t%-8d\t%-8.6g\t%-8.6g\n';
fprintf(1,varNamesSpecs,varNames{:});
outerWCT = tic;
parfor i = 1:I
    runData = data(i,:); % try to communicate only one row to each thread
    % Patches, spans and degrees of freedom:
    if runData{3}.isHybrid && runData{4} > 0
        runData{6} = runData{3}.nonzeroSpanCount;
    else
        runData{6} = 1;
    end
    if isnan(runData{5}) % compute #patches from target #DOFs
        runData{5} = ceil(runData{7}/(runData{4} + runData{6}));
    else
        runData{7} = runData{5}*(runData{4} + runData{6});
    end
    mesh = Mesh(linspace(0,1,runData{5}+1),runData{4},runData{3});
    runData{7} = mesh.dofCount; % update with the actual #DOFs used
    opts = optimset('TolX',runData{9});
    % Assumed characteristic length:
    dl = 1/runData{5}/runData{6}*2.89/(runData{4}^2 + 3*runData{4} + 2.89); % dt = CFL · dl
    tic
    try
        switch runData{8}
            case 'fzero'
                objFun = @(x) dTV_from_dt(runData{3},runData{2},mesh,x);
                [runData{11},runData{12},runData{13}] = fzero(objFun,runData{10}*dl,opts);
            case 'fminbnd'
                objFun = @(x) abs(dTV_from_dt(runData{3},runData{2},mesh,x));
                [runData{11},runData{12},runData{13}] = fminbnd(objFun,runData{10}(1)*dl,runData{10}(2)*dl,opts);
            otherwise
                error('Select either fzero or fminbnd.')
        end
    catch me
        disp(getReport(me,'basic','hyperlinks','on'))
        runData{11} = nan;
        runData{12} = nan;
    end
    runData{14} = toc;
    runData{15} = runData{11}/dl; % CFL = dt / dl
    data(i,:) = runData; % gather back all data
    fprintf(1,varSpecs,runData{1},I,class(runData{2}),class(runData{3}),runData{4:end});
end
disp(' ')
toc(outerWCT)

%% Write to disk
try
    if ~isfile(filename)
        fprintf(1,'Saving as ''%s''.\n',filename);
        fileID = fopen(filename,'w');
    else
        fileID = fopen(filename,'a+');
        fprintf(1,'Appending to ''%s''.\n',filename);
        fprintf(fileID,'\n');
    end
    tic
    fprintf(fileID,varNamesSpecs,varNames{:});
    for i = 1:I
        fprintf(fileID,varSpecs,data{i,1},I,class(data{i,2}),class(data{i,3}),data{i,4:end});
    end
    fclose(fileID);
    toc
catch me
    disp(getReport(me,'basic','hyperlinks','on'))
    fclose all;
end

%% Auxiliary functions
function data = assembleTestTuple(varargin)
    % Assembles all combinations of an arbitrary number of cell arrays of 
    % any number of variables each into a 2D cell array, such that each row 
    % is a combination and each column is an input cell array.
    %
    % Based on the "fundamental principle of counting": 
    % https://stackoverflow.com/questions/3536833/arbitrary-number-of-nested-loops
    %
    % Preallocation:
    N = cellfun('length',varargin); % number of data points in each category
    M = prod(N); % total number of combinations taking one data point from each category
    data = cell(M,1+nargin);
    % Loop over combinations:
    n = ones(size(N)); % set each category's counter
    for i = 1:M
        % Store each combination's counter:
        data{i,1} = i;
        % Loop over categories:
        for j = 1:nargin
            data{i,1+j} = varargin{j}{n(j)}; % select one cell from each category
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

function dTV = dTV_from_dt(method,timeIntegrator,mesh,dt)
    % Increment in total variation as a function of time-step size.
    % Initial condition is hardcoded as a Gaussian pulse.
    method.project(mesh,[],@(x) Functions.gauss(x));
    dTV = mesh.getTotalVariation;
    timeIntegrator.timeNow = 0;
    if timeIntegrator.launchFixedTimeStep(dt,mesh,0,@(t,x) nan)
        dTV = mesh.getTotalVariation - dTV;
    else
        dTV = 1e300;
        %dTV = mesh.getTotalVariation(1) - dTV;
    end
    close
end

function stop = optimplot_cfl(x,optimValues,~)
    stop = false;
    if optimValues.iteration == 0
        hold all
        return
    end
    scatter(optimValues.iteration,x,'filled','b')
    title(['CFL_{max} \approx ', num2str(x,8)])
    xlabel('Iteration')
    ylabel('CFL_{max}')
    setFancyPlotSettings3
end

function stop = optimplot_dtv(~,optimValues,~)
    stop = false;
    if optimValues.iteration == 0
        hold all
        return
    end
    plot(optimValues.iteration-[1 0],[0 0],'--k')
    scatter(optimValues.iteration,optimValues.fval,'filled','r')
    title(['|\DeltaTV| \approx ', num2str(optimValues.fval,8)])
    xlabel('Iteration')
    ylabel('|\DeltaTV|')
    set(gca,'Yscale','log')
    setFancyPlotSettings3
end

function stop = optimplot_residual(x,optimValues,~)
    stop = false;
    if optimValues.iteration == 0
        hold all
        return
    end
    scatter(optimValues.iteration,abs(optimValues.fval),'filled','r')
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    xlabel('Iteration')
    ylabel('Residual ~ \DeltaTV')
    title(['CFL_{max} \approx ', num2str(x,8)])
    setFancyPlotSettings1
end

function stop = optimplot_error(x,optimValues,~)
    stop = false;
    persistent x0;
    if optimValues.iteration == 0
        x0 = optimValues.intervalb;
        hold all
    else
        relErr = abs(1-x0/x);
        scatter(optimValues.iteration,relErr,'filled','b')
        x0 = x;
        title(['CFL_{max} \approx ', num2str(x,8),'  \Rightarrow  dTV \approx ', num2str(optimValues.fval,8)])
    end
    %set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    xlabel('Iteration')
    ylabel('Relative error estimate ~ 1-CFL^*/CFL')
    setFancyPlotSettings1
end

function stop = optimplot_3D(x,optimValues,~)
    stop = false;
    persistent x0 y0 z0;
    if optimValues.iteration == 0
        x0 = 0;
        y0 = nan;
        z0 = nan;
        hold all
        return
    end
    plot3([x0; optimValues.iteration],[y0; x],[z0; optimValues.fval],'-ok','MarkerFaceColor','b')
    x0 = optimValues.iteration;
    y0 = x;
    z0 = optimValues.fval;
    xlabel('Iteration')
    ylabel('Independent variable ~ x')
    zlabel('Objective function ~ f(x)')
    title(['Iteration ' num2str(x0) '; x \approx ' num2str(y0,8) ', f(x) \approx ' num2str(z0,8)])
    setFancyPlotSettings3
    symlog('Z')
    view(70,30)
    drawnow limitrate
end