clc
clear
%close all

% This script computes the modified wavenumber and spatial operator
% eigenvalues (dimensionless and scaled, respectively, by J and dt) of an 
% arbitrary number of discretizations.
% Output is written to a file.

%% Dependencies
addpath('../../Solver')
addpath('../../Basis')

%% Input
inputData = {
% DGSEM:
    'Filename'              'Time'     'Space'         'K'      'beta'     'eigenmodes'     'fig'
%     'mwa_fv_upwind.dat'     SSP_RK3     DGSEM(0)       60       1          inf             1
%     'mwa_fv_center.dat'     SSP_RK3     DGSEM(0)       60       0          inf             1
    'mwa_dgsem_2_all.dat'   SSP_RK3    DGSEM(2)        60       1          inf              nan
    'mwa_dgsem_3_all.dat'   SSP_RK3    DGSEM(3)        60       1          inf              nan
    'mwa_dgsem_2.dat'       SSP_RK3    DGSEM(2)        60       1          1                2
    'mwa_dgsem_3.dat'       SSP_RK3    DGSEM(3)        60       1          1                2
    'mwa_dgsem_4.dat'       SSP_RK3    DGSEM(4)        60       1          1                2
    'mwa_dgsem_5.dat'       SSP_RK3    DGSEM(5)        60       1          1                2
    'mwa_dgsem_6.dat'       SSP_RK3    DGSEM(6)        60       1          1                2
    'mwa_dgsem_7.dat'       SSP_RK3    DGSEM(7)        60       1          1                2
    'mwa_dgsem_8.dat'       SSP_RK3    DGSEM(8)        60       1          1                2
    'mwa_dgsem_9.dat'       SSP_RK3    DGSEM(9)        60       1          1                2
% %     'mwa_dgiga_1.dat'       SSP_RK3     DGIGA(9,2,0)    60       1          inf             4
% %     'mwa_dgiga_2.dat'       SSP_RK3     DGIGA(8,2,0)    60       1          inf             5
% %     'mwa_dgiga_3.dat'       SSP_RK3     DGIGA(7,2,0)    60       1          inf             6
};

%% Setup
tblIn = cell2table(inputData(2:end,2:end),...
    'VariableNames',inputData(1,2:end),...
    'RowNames',inputData(2:end,1));
runData = table2struct(tblIn);
tblIn.Time = arrayfun(@(x) class(x),tblIn.Time,'UniformOutput',false);
tblIn.Space = arrayfun(@(x) x.getName,tblIn.Space,'UniformOutput',false);
disp(tblIn) % pretty print
I = numel(runData);

%% Batch run
for i = 1:I
    fileID = fopen(tblIn(i,:).Row{1},'w');
    try
        % Exact wavenumbers:
        J = runData(i).Space.basisCount; % number of eigenmodes (i.e. basis functions per patch)
        n = -round(runData(i).K*J/2):round(runData(i).K*J/2); % wavemodes (all of them)
        k = 2*pi/runData(i).K*(n); % (nondimensional) wavenumbers
        % Modified wavenumbers:
        kMod = runData(i).Space.getFourierFootprint(runData(i).beta,k);
        kMod = 1i*kMod;
        % Postprocess & export:
        tic, CFL = runData(i).Time.optimizeCFL(runData(i).Space,runData(i).beta);
        fprintf(fileID,'# %s, CFL_max = %.12g (%.6g s)\n',class(runData(i).Time),CFL,toc);
        fprintf(fileID,'# %s; J = %d, %d patches, beta = %f\n',runData(i).Space.getName,runData(i).Space.basisCount,runData(i).K,runData(i).beta);
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','j','n','k_exact','k_real','k_imag','z_real','z_imag');
        if all(isinf(runData(i).eigenmodes))
            runData(i).eigenmodes = 1:J; % store all eigenmodes
        end
        runData(i).eigenmodes(runData(i).eigenmodes > J) = []; % throw away inexistent modes
        for j = runData(i).eigenmodes % loop over requested eigenmodes
            fprintf(fileID,'%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n',[j*ones(size(n)); n; k/J; real(kMod(j,:))/J; imag(kMod(j,:))/J; CFL*imag(kMod(j,:)); CFL*real(kMod(j,:))]);
            fprintf(fileID,'\n');
        end
        fprintf('Run %d of %d completed by worker %d in %.3g s.\n',i,I,get(getCurrentTask(),'ID'),toc)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(getCurrentTask()) || isnan(runData(i).fig)
            continue
        end
        figure(runData(i).fig)
        subplot(2,2,1)
        plot(k/J,real(kMod(runData(i).eigenmodes,:))/J,'.-')
        hold on
        subplot(2,2,3)
        plot(k/J,imag(kMod(runData(i).eigenmodes,:))/J,'.-')
        hold on
        subplot(2,2,[2 4])
        plot(-1i*kMod(runData(i).eigenmodes,:).'*CFL,'.-')
        hold on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch me
       warning('Run %d of %d did not succeed.\n "%s"',i,I,getReport(me))
    end
    fclose(fileID);
end