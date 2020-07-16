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
%   Category    Values                                 Format        Is block?    
    {
    'dt'        num2cell(logspace(-6,-1,25))        '%.12f'       false
    'ErrorL2'   {nan}                               '%.12f'       false
    'Solver'    {'SSP_RK1'}                         '%s'          true
    'L2'        {nan}                               '%.12f'       false
    'TV'        {nan}                               '%.12f'       false
    'Order'     {nan}                               '%.12f'       false
    }
    {
    'dt'        num2cell(logspace(-4,-1,30))        '%.12f'       false
    'ErrorL2'   {nan}                               '%.12f'       false
    'Solver'    {'SSP_RK2'}                         '%s'          true
    'L2'        {nan}                               '%.12f'       false
    'TV'        {nan}                               '%.12f'       false
    'Order'     {nan}                               '%.12f'       false
    }
    {
    'dt'        num2cell(logspace(-3,-1,20))        '%.12f'       false
    'ErrorL2'   {nan}                               '%.12f'       false
    'Solver'    {'SSP_RK3'}                         '%s'          true
    'L2'        {nan}                               '%.12f'       false
    'TV'        {nan}                               '%.12f'       false
    'Order'     {nan}                               '%.12f'       false
    }
    {
    'dt'        num2cell(logspace(-2,-0,25))        '%.12f'       false
    'ErrorL2'   {nan}                               '%.12f'       false
    'Solver'    {'SSP_RK4_5'}                       '%s'          true
    'L2'        {nan}                               '%.12f'       false
    'TV'        {nan}                               '%.12f'       false
    'Order'     {nan}                               '%.12f'       false
    }
    {
    'dt'        num2cell(logspace(-2,-0,25))        '%.12f'       false
    'ErrorL2'   {nan}                               '%.12f'       false
    'Solver'    {'SSP_RK4_10'}                      '%s'          true
    'L2'        {nan}                               '%.12f'       false
    'TV'        {nan}                               '%.12f'       false
    'Order'     {nan}                               '%.12f'       false
    }
    };
fileNames = {
%   Name                        Is active?
    'order_time_RK1.dat'        false
    'order_time_RK2.dat'        false
    'order_time_RK3.dat'        false
    'order_time_RK45.dat'       true
    'order_time_RK410.dat'      true
    };
%exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) 1-sin(pi*x)/(5*pi));
%exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
exactSolution = @(t,x) exp(-9*pi/4*x.^2);
mesh = Mesh(DGIGA(27,3),[-1 1],Periodic(2),1);

%% Loop over files
if numel(inputData) ~= size(fileNames,1)
    error('Inconsistent input.')
else
    I = sum([fileNames{:,2}]);
end
for i = find([fileNames{:,2}])
    %% Distributed batch run
    % Setup parallel write:
    c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);
    spmd
        cpuFile = (fopen(c.Value));
    end
    % Preprocess:
    try
        dataSpecs = [strjoin([{'%d'}; inputData{i}(:,3)],'\t') '\n'];
        fileData = assembleTestTuple(inputData{i}{:,2});
        J = size(fileData,1);
    catch me
        warning('Preprocessing for batch %d of %d failed.\n "%s"',i,I,getReport(me))
    end
    % Run (distributed loop):
    try
        parfor j = 1:J
            tic
            cpuData = fileData(j,:); % distribute data to CPU
            % Solve current run:
            norms = Norm({'ErrorL2','L2','TV'}); % norms to compute
            scheme = str2func(cpuData{4});
            solver = scheme(...
                Advection,[0 2],...
                'timeDelta',cpuData{2},...
                'norm',norms,...
                'exactSolution',exactSolution);
            solver.initialize(mesh)
            solver.launch(mesh)
            % Postprocess:
            [cpuData{[3 5 6]}] = norms.vals;
            fprintf(c.Value,dataSpecs,cpuData{:});
            fprintf('Run %d of %d in batch %d of %d completed by worker %d in %.3g s.\n',j,J,i,I,get(getCurrentTask(),'ID'),toc)
            close gcf
        end
    catch me
        warning('Batch %d of %d crashed.\n "%s"',i,I,getReport(me))
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate the order of convergence:
        tbl.Order = [nan; log10(tbl.ErrorL2(2:end)./tbl.ErrorL2(1:end-1))./log10(tbl.dt(2:end)./tbl.dt(1:end-1))];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fileData = table2cell(tbl(:,2:end));
        J = size(fileData,1);
        dataSpecs = [strjoin(inputData{i}(:,3),'\t') '\n'];
        % Write data into file:
        fileID = fopen(fileNames{i,1},'w');
        try
            fprintf(fileID,'%s\n',strjoin(inputData{i}(:,1),'\t')); % print headers
            for j = 1:J
                fprintf(fileID,dataSpecs,fileData{j,:}); % print one row
                if j < J
                    % Search for possible block changes:
                    for k = find([inputData{i}{:,4}])
                        if contains('%s',inputData{i}{k,3}) % string category
                            if strcmp(fileData{j,k},fileData{j+1,k})
                                continue
                            end
                        else % assume float category
                            if fileData{j,k} == fileData{j+1,k}
                                continue
                            end
                        end
                        fprintf(fileID,'\n'); % start a new block
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