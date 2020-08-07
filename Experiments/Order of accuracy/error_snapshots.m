clc
clear
%close all

% This script solves individual runs one by one (in parallel) and stores 
% the q^h(x) - q(x) function into a data file.
% This function is sampled in an adaptive way (more points are placed in
% regions where gradients are larger).
% Designed to work with Burgers (scalar equation).

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
% DGSEM, batch 5
    'Filename'                              'dt'                    'basis'         'K'
    'error_burgers_hump_dgsem_1.dat'        0.000106081835513945    DGSEM(11)       10 % stable, dt-independent
    'error_burgers_hump_dgsem_2.dat'        0.0151177507061566      DGSEM(11)       10 % marginally stable
    'error_burgers_hump_dgsem_3.dat'        0.0191448197616996      DGSEM(11)       10 % unstable
% DGSEM, batch 9
    'error_burgers_hump_dgsem_4.dat'        1e-4                    DGSEM(4)        21
    'error_burgers_hump_dgsem_5.dat'        1e-4                    DGSEM(4)        25
    'error_burgers_hump_dgsem_6.dat'        1e-4                    DGSEM(4)        30
% DGSEM, batches 10 to 12
    'error_burgers_hump_dgsem_7.dat'        1e-4                    DGSEM(0)        20
    'error_burgers_hump_dgsem_8.dat'        1e-4                    DGSEM(1)        10
    'error_burgers_hump_dgsem_9.dat'        1e-4                    DGSEM(19)       1
    'error_burgers_hump_dgsem_10.dat'       1e-4                    DGSEM(5)        20
    'error_burgers_hump_dgsem_11.dat'       1e-4                    DGSEM(11)       10
    'error_burgers_hump_dgsem_12.dat'       1e-4                    DGSEM(119)      1
% FR, batches 9 to 12
    'error_burgers_hump_fr_1.dat'           1e-4                    FR('min',2)     300
    'error_burgers_hump_fr_2.dat'           1e-4                    FR('Ga',2)      300
    'error_burgers_hump_fr_3.dat'           1e-4                    FR('LumpLo',2)  300
    'error_burgers_hump_fr_4.dat'           1e-4                    FR('max',2)     300
% FR, batch 17
    'error_burgers_hump_fr_5.dat'           1e-4                    FR('min',4)     74
    'error_burgers_hump_fr_6.dat'           1e-4                    FR('min',4)     88
    'error_burgers_hump_fr_7.dat'           1e-4                    FR('min',4)     105
    'error_burgers_hump_fr_8.dat'           1e-4                    FR('min',4)     126
% DGIGA, batches 46 to 49 
    'error_burgers_hump_dgiga_1.dat'        1e-4                    DGIGA(2,2,1)    15
    'error_burgers_hump_dgiga_2.dat'        1e-4                    DGIGA(2,3,1)    10
    'error_burgers_hump_dgiga_3.dat'        1e-4                    DGIGA(4,2,1)    10
    'error_burgers_hump_dgiga_4.dat'        1e-4                    DGIGA(2,4,3)    10
% DGIGA, batches 47 and 49
    'error_burgers_hump_dgiga_5.dat'        1e-4                    DGIGA(2,5,1)    10
    'error_burgers_hump_dgiga_6.dat'        1e-4                    DGIGA(2,6,1)    10
    'error_burgers_hump_dgiga_7.dat'        1e-4                    DGIGA(2,8,7)    10
    'error_burgers_hump_dgiga_8.dat'        1e-4                    DGIGA(2,10,9)   10
% DGIGA, batch 70
    'error_burgers_hump_dgiga_9.dat'        1e-4                    DGIGA(16,4,3)   7 % unstable
    'error_burgers_hump_dgiga_10.dat'       1e-4                    DGIGA(16,4,3)   11 % unstable
    'error_burgers_hump_dgiga_11.dat'       1e-4                    DGIGA(16,4,3)   45 % stable?
% DGIGA, batches 79 and 80
    'error_burgers_hump_dgiga_12.dat'       1e-4                    DGIGA(53,2,1)   1 % stable
    'error_burgers_hump_dgiga_13.dat'       1e-4                    DGIGA(53,3,2)   1 % unstable
% DGIGA, batch 106
    'error_burgers_hump_dgiga_14.dat'       1e-4                    DGIGA(8,15,0)   1 % stable?
    'error_burgers_hump_dgiga_15.dat'       1e-4                    DGIGA(8,16,0)   1 % unstable
% DGIGA, batch 114
    'error_burgers_hump_dgiga_16.dat'       1e-4                    DGIGA(8,5,4)   20 % stable
    'error_burgers_hump_dgiga_17.dat'       1e-4                    DGIGA(8,6,5)   20 % unstable
};
exactSolution = @(t,x) smoothBurgersExact(t,x,@(x) exp(-9*pi/4*x.^2));
N = 1801; % total number of (distinct) sample locations

%% Setup
tblIn = cell2table(inputData(2:end,2:end),...
    'VariableNames',inputData(1,2:end),...
    'RowNames',inputData(2:end,1));
runData = table2struct(tblIn);
if any(any(~cellfun(@isscalar,inputData(2:end,2:end)))) % safety check
    error('Please, specify one run per file only.')
end
tblIn.basis = arrayfun(@(x) x.getName,tblIn.basis,'UniformOutput',false);
disp(tblIn) % pretty print
J = numel(runData);

%% Parallel run
parfor j = 1:J
    fileID = fopen(tblIn(j,:).Row{1},'w');
    try
        norms = Norm({'TV','L2','ErrorL2'}); % norms to compute
        mesh = Mesh(runData(j).basis,[-1 1],Periodic(2),runData(j).K);
        solver = SSP_RK3(...
            Burgers,[0 .25],...
            'timeDelta',runData(j).dt,...
            'norm',norms,...
            'exactSolution',exactSolution);
        solver.initialize(mesh)
        solver.launch(mesh)
        close gcf
        if solver.timeNow ~= solver.timeStop % run was not succesful
            fprintf('Run %d of %d has been skipped (worker %d, %.3g s).\n',j,J,get(getCurrentTask(),'ID'),toc)
            continue % move on to the next run
        end
        % Postprocess & export:
        fprintf(fileID,'# dt = %.12g\n',runData(j).dt);
        fprintf(fileID,'# %s; %d elements, N_dofs = %d\n',mesh.bases.getName,runData(j).K,mesh.dofCount);
        fprintf(fileID,'# TV = %.12g, L2 = %.12g, ErrorL2 = %.12g\n',norms.vals);
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\n','x','error','exact','approx','residual','initial_exact');
        n = round((N-1)/mesh.elementCount+1); % number of samples per element
        for element = mesh.elements
            x = linspace(element.edgeL.coord,element.edgeR.coord,n);
            z = [
                exactSolution(solver.timeNow,x)
                element.interpolateStateAtCoords(x)
            	element.interpolateResidualAtCoords(x)
                exactSolution(solver.timeStart,x)
            ];
            fprintf(fileID,'%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n',[x; abs(z(1,:)-z(2,:)); z]);
            fprintf(fileID,'\n'); % blank line
        end
        fprintf('Run %d of %d completed by worker %d in %.3g s.\n',j,J,get(getCurrentTask(),'ID'),toc)
    catch me
       warning('Run %d of %d did not succeed.\n "%s"',j,J,getReport(me))
    end
    fclose(fileID);
end