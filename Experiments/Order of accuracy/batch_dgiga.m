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

%% Input
maxDOFs = 900; % will skip any runs beyond this number
inputData = {
    'Filename'                  'dt'        'Method'                    'K'                             'p'   'k'   'kappa'
    'order_dgiga_k_1.dat'       1e-4        {'DGIGA','DGIGA_nodal'}     logspacei(1,300,25)             2     1     0
    'order_dgiga_k_2.dat'       1e-4        {'DGIGA','DGIGA_nodal'}     logspacei(1,225,25)             3     1     0
    'order_dgiga_k_3.dat'       1e-4        {'DGIGA','DGIGA_nodal'}     logspacei(1,180,25)             4     1     0
    'order_dgiga_k_4.dat'       1e-4        'DGIGA'                     logspacei(1,150,25)             2     4     0
    'order_dgiga_k_5.dat'       1e-4        'DGIGA'                     logspacei(1,100,25)             2     4     1
    'order_dgiga_k_6.dat'       1e-4        'DGIGA'                     logspacei(1,27,25)              2     16    0
    'order_dgiga_k_7.dat'       1e-4        'DGIGA'                     logspacei(1,50,25)              2     16    1
    'order_dgiga_k_8.dat'       1e-4        'DGIGA'                     logspacei(1,69,25)              3     4     0
    'order_dgiga_k_9.dat'       1e-4        'DGIGA'                     logspacei(1,128,25)             3     4     2
    'order_dgiga_k_10.dat'      1e-4        'DGIGA'                     1:18                            3     16    0
    'order_dgiga_k_11.dat'      1e-4        'DGIGA'                     logspacei(1,47,25)              3     16    2
    'order_dgiga_k_12.dat'      1e-4        'DGIGA'                     logspacei(1,52,25)              4     4     0
    'order_dgiga_k_13.dat'      1e-4        'DGIGA'                     logspacei(1,112,25)             4     4     3
    'order_dgiga_k_14.dat'      1e-4        'DGIGA'                     1:13                            4     16    0
    'order_dgiga_k_15.dat'      1e-4        'DGIGA'                     logspacei(1,45,25)              4     16    3
};
%exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) 1-sin(pi*x)*2/(5*pi));
 exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
%exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) cos(.5*pi*x).^5);
%exactSolution = @(t,x) exp(-9*pi/4*x.^2);

%% Setup
inputTable = cell2table(inputData(2:end,2:end),...
    'VariableNames',inputData(1,2:end),...
    'RowNames',inputData(2:end,1)); % pretty print
disp(inputTable)
I = size(inputTable,1);

%% Loop over batches
for i = 1:I
    %% Distributed batch run
    % Setup parallel write:
    c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);
    spmd
        cpuFile = (fopen(c.Value));
        fprintf(c.Value,'Run\tdt\tMethod\tDOFs\tK\tp\tk\tkappa\tWCT\tTV\tL2\tErrorL2\n');
    end
    % Preprocess:
    try
        for k = size(inputTable,2):-1:1 % column-wise num2cell
            if iscell(inputTable{i,k})
                if ischar(inputTable{i,k}{1})
                    batchData{k} = inputTable{i,k};
                    continue
                end
                batchData{k} = num2cell(inputTable{i,k}{1});
            else
                batchData{k} = num2cell(inputTable{i,k});
            end
        end
        runData = cell2struct(assembleTestTuple(batchData{:}),['Run' inputTable.Properties.VariableNames],2); % combines batch into individual runs
        %%%
        % Replace Inf by p-1 (smoothness):
        aux = num2cell([runData.p]-1);
        [runData(isinf([runData.kappa]) | isnan([runData.kappa])).kappa] = aux{:};
        % Filter out "incompatible" runs:
        runData([runData.kappa] >= [runData.p],:) = []; % smoothness > degree-1, skip!
        % Filter out "excessive" runs:
        runData([runData.K].*([runData.p].*[runData.k] - [runData.kappa].*[runData.k] + [runData.kappa] + 1) > maxDOFs,:) = []; % too many DOFs, skip!
        %%%
        J = numel(runData);
    catch me
        warning('Preprocessing for batch %d of %d failed.\n "%s"',i,I,getReport(me))
    end
    % Run (distributed loop):
    try
        parfor j = 1:J
            tic
            % Solve current run:
            norms = Norm({'TV','L2','ErrorL2'}); % norms to compute
            method = str2func(runData(j).Method);
            mesh = Mesh(method(runData(j).k,runData(j).p,runData(j).kappa),[-1 1],Periodic(2),runData(j).K);
            solver = SSP_RK3(...
                Burgers,[0 .25],...
                'timeDelta',runData(j).dt,...
                'norm',norms,...
                'exactSolution',exactSolution);
            solver.initialize(mesh)
            solver.launch(mesh)
            % Postprocess:
            close gcf
            if solver.timeNow ~= solver.timeStop % run was not succesful
                fprintf('Run %d of %d in batch %d of %d has been skipped (worker %d, %.3g s).\n',j,J,i,I,get(getCurrentTask(),'ID'),toc)
                continue % do not write it
            end
            outputData = {runData(j).Run runData(j).dt runData(j).Method mesh.dofCount runData(j).K runData(j).p runData(j).k mesh.bases.smoothness solver.wallClockTime norms(:).vals};
            fprintf(c.Value,'%d\t%.12f\t%s\t%d\t%d\t%d\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\n',outputData{:});
            fprintf('Run %d of %d in batch %d of %d completed by worker %d in %.3g s.\n',j,J,i,I,get(getCurrentTask(),'ID'),toc)
        end
    catch me
        warning('Batch %d of %d crashed.\n "%s"',i,I,getReport(me))
    end
    clear c % free memory and close temporary files
    %% Export to file
    try
        % Read temporary files:
        spmd
            cpuTable = readtable(cpuFile,'Delimiter','\t','ReadVariableNames',true);
        end
        outputTable = sortrows(vertcat(cpuTable{:}));
        clear cpuTable % free memory
        % Write to output file:
        try
            writetable(outputTable(:,2:end),inputTable.Row{i},'Delimiter','\t')
            spmd
                delete(cpuFile) % delete them
            end
        catch me
            warning('Export to file for batch %d of %d failed.\n "%s"',i,I,getReport(me))
        end
    catch me
        warning('Could not read temporary files for batch %d of %d.\n "%s"',i,I,getReport(me))
    end
end