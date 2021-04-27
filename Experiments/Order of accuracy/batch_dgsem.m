clc
clear

%% Input
inputData = {
    'Filename'              'dt'                'K'                 'p' 
    'order_dgsem_dt_1.dat'  logspace(-5,-2,30)  300                 2
    'order_dgsem_dt_2.dat'  logspace(-5,-2,30)  225                 3
    'order_dgsem_dt_3.dat'  logspace(-5,-2,30)  180                 4
    'order_dgsem_dt_4.dat'  logspace(-5,-1,40)  20                  5
    'order_dgsem_dt_5.dat'  logspace(-5,-1,40)  10                  11
    'order_dgsem_dt_6.dat'  logspace(-5,-1,40)  1                   119
    'order_dgsem_k_1.dat'   1e-4                logspacei(1,300,25) 2
    'order_dgsem_k_2.dat'   1e-4                logspacei(1,225,25) 3
    'order_dgsem_k_3.dat'   1e-4                logspacei(1,180,25) 4
    'order_dgsem_p_1.dat'   1e-4                20                  0:5
    'order_dgsem_p_2.dat'   1e-4                10                  0:11
    'order_dgsem_p_3.dat'   1e-4                1                   [0 logspacei(1,119,24)]
};
exactSolution = @(t,x) Burgers.MOC(t,x,@(x) exp(-9*pi/4*x.^2),[-1 1]);

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
        fprintf(c.Value,'Run\tdt\tDOFs\tK\tp\tWCT\tTV\tL2\tErrorL2\n');
    end
    % Preprocess:
    try
        for k = size(inputTable,2):-1:1 % column-wise num2cell
            if iscell(inputTable{i,k})
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
            mesh = Mesh(DGSEM(runData(j).p),[-1 1],Periodic(2),runData(j).K);
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
            outputData = {runData(j).Run runData(j).dt mesh.dofCount runData(j).K runData(j).p solver.wallClockTime norms(:).vals};
            fprintf(c.Value,'%d\t%.12f\t%d\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\n',outputData{:});
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