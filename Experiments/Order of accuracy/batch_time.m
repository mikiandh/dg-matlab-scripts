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
%   Category    Values                          Format          Is block?
    'CFL'       num2cell(logspace(-3,-1,5))     '%.6f'          false
    'ErrorL2'   {nan}                           '%.6f'          false
    'Solver'    {'SSP_RK1','SSP_RK2','SSP_RK3'} '%s'            true
    'L2'        {nan}                           '%.6f'          false
    'TV'        {nan}                           '%.6f'          false
    }
    {
    %   Category    Values                          Format          Is block?
    'CFL'       num2cell(logspace(-2,-0,10))    '%.6f'          false
    'ErrorL2'   {nan}                           '%.6f'          false
    'Solver'    {'SSP_RK4_10'}                  '%s'            true
    'L2'        {nan}                           '%.6f'          false
    'TV'        {nan}                           '%.6f'          false    
    }
    };
fileNames = {
    'order_time_DGSEM_1.dat';
    'order_time_DGSEM_2.dat';
    };
exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
mesh = Mesh(DGSEM(2),[-1 1],Periodic(2),5);

%% Loop over files
if numel(inputData) ~= numel(fileNames)
    error('Inconsistent input.')
else
    I = numel(inputData);
end
for i = 1:I
    %% Distributed batch run
    % Setup parallel write:
    c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);
    spmd
        cpuFile = (fopen(c.Value));
    end
    % Preprocess:
    dataSpecs = ['%d,' strjoin(inputData{i}(:,3),',') '\n'];
    fileData = assembleTestTuple(inputData{i}{:,2});
    J = size(fileData,1);
    % Run (distributed loop):
    try
        parfor j = 1:J
            tic
            cpuData = fileData(j,:); % distribute data to CPU
            % Solve current run:
            norms = Norm({'ErrorL2','L2','TV'}); % norms to compute
            scheme = str2func(cpuData{4});
            solver = scheme(...
                Burgers,[0 .3],...
                'courantNumber',cpuData{2},...
                'norm',norms,...
                'exactSolution',exactSolution);
            solver.initialize(mesh)
            solver.launch(mesh)
            % Postprocess:
            [cpuData{[3 5 6]}] = norms.vals;
            fprintf(c.Value,dataSpecs,cpuData{:});
            fprintf('Worker %d: run %d of %d in batch %d of %d completed (%.3g s).\n',get(getCurrentTask(),'ID'),j,J,i,I,toc)
            close gcf
            %%%fileData(j,:) = cpuData;
        end
    catch me
        warning('Batch %d of %d crashed:\n "%s"',i,I,getReport(me))
    end
    clear c % close temporary files
    %% Export to file
    try
        % Read temporary files:
        spmd
            tbl0 = readtable(cpuFile,'ReadVariableNames',false);
        end
        tbl = sortrows(vertcat(tbl0{:}));
        tbl.Properties.VariableNames = [{'Run'}; inputData{i}(:,1)];
        fileData = table2cell(tbl(:,2:end));
        % Write data into file:
        fileID = fopen(fileNames{i},'w');
        try
            fprintf(fileID,'%s\n',strjoin(inputData{i}(:,1),',')); % print headers
            for j = 1:J
                fprintf(fileID,dataSpecs,fileData{j,:}); % print one row
                if j < J
                    % Search for possible block changes:
                    for k = find([inputData{i}{:,4}])
                        if any(fileData{j,k} ~= fileData{j+1,k})
                            fprintf(fileID,'\n'); % start a new block
                        end
                    end
                end
            end
            spmd
                delete(cpuFile) % delete them
            end
        catch me
            warning('Export to file for batch %d of %d failed.\n "%s"',i,I,getReport(me))
        end
        fclose(fileID);
    catch me
        warning('Could not read temporary files for batch %d of %d.\n "%s"',i,I,getReport(me))
    end
end