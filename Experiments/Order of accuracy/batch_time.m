clc
clear
%close all

%% Dependencies
addpath('../../Limiting')
addpath('../../Physics')
addpath('../../Solver')
addpath('../../Basis')
addpath('../../Mesh')
addpath('../../Math')
addpath('../../Extra')

%% Setup
inputData = {
%   Time convergence
    {
%   Category    Values                      Format          Is block?
    'CFL'       num2cell(logspace(-3,-1,5)) '%.6g'          false
    'ErrorL2'   {nan}                       '%.6g'          false
    'Solver'    {'SSP_RK1','SSP_RK2'}       '%s'            true
    'L2'        {nan}                       '%.6g'          false
    'TV'        {nan}                       '%.6g'          false
    }
    };
fileNames = {
    'order_time_DGSEM.dat';
    };
exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
mesh = Mesh(DGSEM(2),[-1 1],Periodic(2),5);

%% Distributed batch run
% Consistency check:
if numel(inputData) ~= numel(fileNames)
    error('Inconsistent input.')
else
    I = numel(inputData);
end
% Loop over files:
for i = 1:I
    % Preprocess:
    fileData = assembleTestTupleBlocks(inputData{i}{:,2});
    J = size(fileData,1);
    % Run (distributed loop):
    for j = 1:J
        tic
        cpuData = fileData(j,3:end); % distribute data to CPU
        % Solve current run:
        norms = Norm({'ErrorL2','L2','TV'}); % norms to compute
        scheme = str2func(cpuData{3});
        solver = scheme(...
            Burgers,[0 .3],...
            'courantNumber',cpuData{1},...
            'norm',norms,...
            'exactSolution',exactSolution);
        solver.initialize(mesh)
        solver.launch(mesh)
        % Postprocess:
        [cpuData{[2 4 5]}] = norms.vals;
        close gcf
        fileData(j,3:end) = cpuData; % gather data back
        % Print status:
        fprintf('Worker %d, run %d of %d done (%.3g s); block %d of %d.\n',get(getCurrentTask(),'ID'),j,J,toc,i,I)
    end
    % Export to file:
    fileID = fopen(fileNames{i},'w');
    try
    fprintf(fileID,'Block,Run,%s\n',strjoin(inputData{i}(:,1),',')); % header format
    dataSpecs = ['%d,%d,' strjoin(inputData{i}(:,3),',') '\n']; % data format
    for j = 1:J
        fprintf(fileID,dataSpecs,fileData{j,:}); % print one row
        if j < J
            m = fileData{j+1,1};
            if m ~= fileData{j,1} && inputData{i}{m,4}
                fprintf(fileID,'\n'); % start a new block
            end
        end
    end
    catch
        warning('Export to file failed.')
    end
    fclose(fileID);
end