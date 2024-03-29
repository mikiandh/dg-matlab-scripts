clc
clear

%% Input
maxDOFs = 950; % will skip any runs beyond this number
inputData = {
    'Filename'                  'dt'                'Method'                    'K'                   'p'     'k'                      'kappa'
% DG h-refinement -------------------------------------------------------------------------------------------------------------------------------------
    'order_dgiga_dt_1.dat'      logspace(-5,-2,30)  'DGIGA'                     300            		  2       1   					 	1 
    'order_dgiga_dt_2.dat'      logspace(-5,-2,30)  'DGIGA'                     225            		  3       1                      	2 
    'order_dgiga_dt_3.dat'      logspace(-5,-2,30)  'DGIGA'                 	  180            		  4       1                      	3 
    'order_dgiga_dt_4.dat'      logspace(-5,-2,30)  'DGIGA'                     100             	    2  	   4                         0 
    'order_dgiga_dt_5.dat'      logspace(-5,-2,30)  'DGIGA'                     150                     2      4                         1 
    'order_dgiga_dt_6.dat'      logspace(-5,-2,30)  'DGIGA'                     27				 	  2   	  16                        0
    'order_dgiga_dt_7.dat'      logspace(-5,-2,30)  'DGIGA'                     50				 	  2   	  16                        1
    'order_dgiga_dt_8.dat'      logspace(-5,-2,30)  'DGIGA'                     69				 	  3   	  4                         0
    'order_dgiga_dt_9.dat'      logspace(-5,-2,30)  'DGIGA'                     128					  3   	  4                         2
    'order_dgiga_dt_10.dat'     logspace(-5,-2,30)  'DGIGA'                     18               	      3   	  16                        0
    'order_dgiga_dt_11.dat'     logspace(-5,-2,30)  'DGIGA'                     47 					  3   	  16                        2
    'order_dgiga_dt_12.dat'     logspace(-5,-2,30)  'DGIGA'                     52 	  				  4   	  4                         0
    'order_dgiga_dt_13.dat'     logspace(-5,-2,30)  'DGIGA'                     112	  				  4   	  4                         3
    'order_dgiga_dt_14.dat'     logspace(-5,-2,30)  'DGIGA'                     13                      4   	  16                        0
    'order_dgiga_dt_15.dat'     logspace(-5,-2,30)  'DGIGA'                     45				 	  4   	  16                        3
%     
    'order_dgiga_k_1.dat'       1e-4                {'DGIGA','DGIGA_nodal'}     logspacei(1,300,25)	  2   	  1                         0
    'order_dgiga_k_2.dat'       1e-4                {'DGIGA','DGIGA_nodal'}     logspacei(1,225,25)	  3   	  1                         0
    'order_dgiga_k_3.dat'       1e-4                {'DGIGA','DGIGA_nodal'}     logspacei(1,180,25)	  4   	  1                         0
    'order_dgiga_k_4.dat'       1e-4                'DGIGA'                     logspacei(1,100,25)	  2   	  4                         0
    'order_dgiga_k_5.dat'       1e-4                'DGIGA'                     logspacei(1,150,25)	  2   	  4                         1
    'order_dgiga_k_6.dat'       1e-4                'DGIGA'                     logspacei(1,27,25) 	  2   	  16                        0
    'order_dgiga_k_7.dat'       1e-4                'DGIGA'                     logspacei(1,50,25) 	  2   	  16                        1
    'order_dgiga_k_8.dat'       1e-4                'DGIGA'                     logspacei(1,69,25) 	  3   	  4                         0
    'order_dgiga_k_9.dat'       1e-4                'DGIGA'                     logspacei(1,128,25)	  3   	  4                         2
    'order_dgiga_k_10.dat'      1e-4                'DGIGA'                     1:18               	  3   	  16                        0
    'order_dgiga_k_11.dat'      1e-4                'DGIGA'                     logspacei(1,47,25) 	  3   	  16                        2
    'order_dgiga_k_12.dat'      1e-4                'DGIGA'                     logspacei(1,52,25) 	  4   	  4                         0
    'order_dgiga_k_13.dat'      1e-4                'DGIGA'                     logspacei(1,112,25)	  4   	  4                         3
    'order_dgiga_k_14.dat'      1e-4                'DGIGA'                     1:13               	  4   	  16                        0
    'order_dgiga_k_15.dat'      1e-4                'DGIGA'                     logspacei(1,45,25) 	  4   	  16                        3
% p-refinement ----------------------------------------------------------------------------------------------------------------------------------------
	'order_dgiga_dt_16.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  16      1   					 	0
	'order_dgiga_dt_17.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  15      2   					 	0
	'order_dgiga_dt_18.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  11      4   					 	0
	'order_dgiga_dt_19.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  5       8   					 	0
	'order_dgiga_dt_20.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  12      2   					 	inf
	'order_dgiga_dt_21.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  11      4   					 	inf
	'order_dgiga_dt_22.dat'     logspace(-5,-2,30)  'DGIGA'                     20            		  6       8   					 	inf
%	
    'order_dgiga_p_1.dat'       1e-4                'DGIGA'                     [1 20]             	  1:49	  2.^(0:3)                  0
    'order_dgiga_p_2.dat'       1e-4                'DGIGA'                     [10 20]            	  1:12	  2.^(0:3)                  inf
% DGIGA knot insertion --------------------------------------------------------------------------------------------------------------------------------
	'order_dgiga_dt_23.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  2   	  449       				0
    'order_dgiga_dt_24.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  2   	  898				        inf
    'order_dgiga_dt_25.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  2   	  44				        0
    'order_dgiga_dt_26.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  2   	  87				        inf
    'order_dgiga_dt_27.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  2   	  22	                    0
    'order_dgiga_dt_28.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  2   	  43				        inf
    'order_dgiga_dt_29.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  3   	  299				        0
    'order_dgiga_dt_30.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  3   	  53                        inf
    'order_dgiga_dt_31.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  3   	  26                        0
    'order_dgiga_dt_32.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  3   	  33                        inf
    'order_dgiga_dt_33.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  3   	  14                        0
    'order_dgiga_dt_34.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  3   	  35                        inf
    'order_dgiga_dt_35.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  4   	  224                       0
    'order_dgiga_dt_36.dat'     logspace(-5,-2,30)  'DGIGA'                     1                  	  4   	  9                         inf
    'order_dgiga_dt_37.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  4   	  22    	                0
    'order_dgiga_dt_38.dat'     logspace(-5,-2,30)  'DGIGA'                     10                 	  4   	  16        				inf
    'order_dgiga_dt_39.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  4   	  11                      	0
    'order_dgiga_dt_40.dat'     logspace(-5,-2,30)  'DGIGA'                     20                 	  4   	  16        				inf
%
    'order_dgiga_knot_1.dat'    1e-4                'DGIGA'                     1                  	  2   	  logspacei(1,449,25)       0
    'order_dgiga_knot_2.dat'    1e-4                'DGIGA'                     1                  	  2   	  logspacei(1,898,25)       1
    'order_dgiga_knot_3.dat'    1e-4                'DGIGA'                     10                 	  2   	  logspacei(1,44,25)       0
    'order_dgiga_knot_4.dat'    1e-4                'DGIGA'                     10                 	  2   	  logspacei(1,87,25)        1
    'order_dgiga_knot_5.dat'    1e-4                'DGIGA'                     20                 	  2   	  1:22                      0
    'order_dgiga_knot_6.dat'    1e-4                'DGIGA'                     20                 	  2   	  logspacei(1,43,25)        1
    'order_dgiga_knot_7.dat'    1e-4                'DGIGA'                     1                  	  3   	  logspacei(1,299,25)       0
    'order_dgiga_knot_8.dat'    1e-4                'DGIGA'                     1                  	  3   	  logspacei(1,897,25)       2
    'order_dgiga_knot_9.dat'    1e-4                'DGIGA'                     10                 	  3   	  logspacei(1,299,25)       0
    'order_dgiga_knot_10.dat'   1e-4                'DGIGA'                     10                 	  3   	  logspacei(1,87,25)        2
    'order_dgiga_knot_11.dat'   1e-4                'DGIGA'                     20                 	  3   	  logspacei(1,29,25)        0
    'order_dgiga_knot_12.dat'   1e-4                'DGIGA'                     20                 	  3   	  logspacei(1,42,25)        2
    'order_dgiga_knot_13.dat'   1e-4                'DGIGA'                     1                  	  4   	  logspacei(1,224,25)       0
    'order_dgiga_knot_14.dat'   1e-4                'DGIGA'                     1                  	  4   	  logspacei(1,896,25)       3
    'order_dgiga_knot_15.dat'   1e-4                'DGIGA'                     10                 	  4   	  1:22                      0
    'order_dgiga_knot_16.dat'   1e-4                'DGIGA'                     10                 	  4   	  logspacei(1,86,25)        3
    'order_dgiga_knot_17.dat'   1e-4                'DGIGA'                     20                 	  4   	  1:11                      0
    'order_dgiga_knot_18.dat'   1e-4                'DGIGA'                     20                 	  4   	  logspacei(1,41,25)        3
% 4 refinement directions ---------------------------------------------------------------------------------------------------------------------------- 
    'order_dgiga_dt_41.dat'     logspace(-5,-2,30)  'DGIGA'                     10                      2       2                         1
    'order_dgiga_dt_42.dat'     logspace(-5,-2,30)  'DGIGA'                     30                      2       2                         1
    'order_dgiga_dt_43.dat'     logspace(-5,-2,30)  'DGIGA'                     10                      6       2                         1
    'order_dgiga_dt_44.dat'     logspace(-5,-2,30)  'DGIGA'                     10                      2       10                        1
    'order_dgiga_dt_45.dat'     logspace(-5,-2,30)  'DGIGA'                     10                      10      2                         inf
%
    'order_dgiga_extra_1.dat'   1e-4                'DGIGA'                     10:30                   2       2                         1
    'order_dgiga_extra_2.dat'   1e-4                'DGIGA'                     10                      2:6     2                         1
    'order_dgiga_extra_3.dat'   1e-4                'DGIGA'                     10                      2       2:10                      1
    'order_dgiga_extra_4.dat'   1e-4                'DGIGA'                     10                      2:10    2                         inf
% IGA & CG knot insertion -----------------------------------------------------------------------------------------------------------------------------
    'order_iga_knot_1.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       1       logspacei(1,899,25)       0
    'order_iga_knot_2.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       2       logspacei(1,449,25)       0
    'order_iga_knot_3.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       2       logspacei(1,898,25)       1
    'order_iga_knot_4.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       3       logspacei(1,299,25)       0
    'order_iga_knot_5.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       3       logspacei(1,897,25)       2
    'order_iga_knot_6.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       4       logspacei(1,224,25)       0
    'order_iga_knot_7.dat'      1e-4                {'DGIGA','DGIGA_nodal'}     1                       4       logspacei(1,896,25)       3
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