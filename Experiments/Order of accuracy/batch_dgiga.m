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
inputData = {
    'Filename'               'dt'                'K'            'p'      'k'    'kappa'
% K <-> p:
    'order_dgiga_dt_1.dat'   logspace(-5,-2,30)  300            2        1      1 % DOF = 900
    'order_dgiga_dt_2.dat'   logspace(-5,-2,30)  225            3        1      2 % DOF = 900
    'order_dgiga_dt_3.dat'   logspace(-5,-2,30)  180            4        1      3 % DOF = 900
% p <-> k:
    'order_dgiga_dt_4.dat'   logspace(-5,-2,30)  10             2        10     1 % DOF = 120
    'order_dgiga_dt_5.dat'   logspace(-5,-2,30)  10             3        9      2 % DOF = 120
    'order_dgiga_dt_6.dat'   logspace(-5,-2,30)  10             4        8      3 % DOF = 120
    'order_dgiga_dt_7.dat'   logspace(-5,-2,30)  20             2        4      1 % DOF = 120
    'order_dgiga_dt_8.dat'   logspace(-5,-2,30)  20             3        3      2 % DOF = 120
    'order_dgiga_dt_9.dat'   logspace(-5,-2,30)  20             4        2      3 % DOF = 120
% k <-> kappa:
    'order_dgiga_dt_10.dat'  logspace(-5,-2,30)  25             3        33     2 % DOF = 900
    'order_dgiga_dt_11.dat'  logspace(-5,-2,30)  25             3        17     1 % DOF = 900
    'order_dgiga_dt_12.dat'  logspace(-5,-2,30)  25             3        12     0 % DOF = 925
% K <-> k:
    'order_dgiga_dt_13.dat'  logspace(-5,-2,30)  75             2        10     1 % DOF = 900
    'order_dgiga_dt_14.dat'  logspace(-5,-2,30)  36             2        23     1 % DOF = 900
    'order_dgiga_dt_15.dat'  logspace(-5,-2,30)  20             2        43     1 % DOF = 900
% K <-> kappa:
    'order_dgiga_dt_16.dat'  logspace(-5,-2,30)  30             3        27     2 % DOF = 900
    'order_dgiga_dt_17.dat'  logspace(-5,-2,30)  16             3        27     1 % DOF = 896
    'order_dgiga_dt_18.dat'  logspace(-5,-2,30)  11             3        27     0 % DOF = 902
};
%exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) 1-sin(pi*x)*2/(5*pi));
exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
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
        fprintf(c.Value,'Run\tdt\tDOFs\tK\tp\tk\tkappa\tWCT\tTV\tL2\tErrorL2\n');
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
            mesh = Mesh(DGIGA(runData(j).k,runData(j).p,runData(j).kappa),[-1 1],Periodic(2),runData(j).K);
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
            outputData = {runData(j).Run runData(j).dt mesh.dofCount runData(j).K runData(j).p runData(j).k mesh.bases.smoothness solver.wallClockTime norms(:).vals};
            fprintf(c.Value,'%d\t%.12f\t%d\t%d\t%d\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\n',outputData{:});
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
        % Estimate the order of convergence (in space):
        outputTable.Order = [nan; log10(outputTable.ErrorL2(2:end)./outputTable.ErrorL2(1:end-1))./log10(outputTable.DOFs(1:end-1)./outputTable.DOFs(2:end))];
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