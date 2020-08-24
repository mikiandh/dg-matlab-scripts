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
    'Filename'              'Time'     'Space'         'K'      'beta'     'eigenmodes'    'fig'
%     'mwa_dgsem_0.dat'       SSP_RK3    DGSEM(0)        60       1          inf             nan
%     'mwa_dgsem_1.dat'       SSP_RK3    DGSEM(1)        60       1          inf             nan
%     'mwa_dgsem_2.dat'       SSP_RK3    DGSEM(2)        60       1          inf             nan
%     'mwa_dgsem_3.dat'       SSP_RK3    DGSEM(3)        60       1          inf             nan
%     'mwa_dgsem_4.dat'       SSP_RK3    DGSEM(4)        60       1          inf             nan
%     'mwa_dgsem_5.dat'       SSP_RK3    DGSEM(5)        60       1          inf             nan
%     'mwa_dgsem_6.dat'       SSP_RK3    DGSEM(6)        60       1          inf             nan
%     'mwa_dgsem_7.dat'       SSP_RK3    DGSEM(7)        60       1          inf             nan
%     'mwa_dgsem_8.dat'       SSP_RK3    DGSEM(8)        60       1          inf             nan
%     'mwa_dgsem_9.dat'       SSP_RK3    DGSEM(9)        60       1          inf             nan
%     'mwa_dgsem_21.dat'      SSP_RK3    DGSEM(21)       60       1          inf             nan
%     'mwa_dgsem_50.dat'      SSP_RK3    DGSEM(50)       60       1          inf             nan
%     'mwa_dgsem_119.dat'     SSP_RK3    DGSEM(119)      60       1          inf             nan
    'mwa_dgiga_1.dat'       SSP_RK4_10    DGIGA(9,2,0)    60       1          inf             4
    'mwa_dgiga_2.dat'       SSP_RK4_10    DGIGA(8,2,0)    60       1          inf             5
    'mwa_dgiga_3.dat'       SSP_RK4_10    DGIGA(7,2,0)    60       1          inf             6
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
        % Wavenumbers and stuff:
        J = runData(i).Space.basisCount; % number of eigenmodes (i.e. basis functions per patch)
        N = round(runData(i).K*J/2); % highest wavemode
        k = 2*pi/runData(i).K*(-N:N); % all (nondimensional) wavenumbers
        kMod = nan(J,2*N+1); % corresponding modified wavenumbers
        % Fourier residual matrices:
        if isa(runData(i).Space,'FR')
            C = -runData(i).Space.gradientMatrix';
            E = (1+runData(i).beta)*runData(i).Space.correctionsL'*runData(i).Space.left' + (1-runData(i).beta)*runData(i).Space.correctionsR'*runData(i).Space.right';
            Eneg = -(1+runData(i).beta)*runData(i).Space.correctionsL'*runData(i).Space.right';
            Epos = -(1-runData(i).beta)*runData(i).Space.correctionsR'*runData(i).Space.left';
        else
            C = runData(i).Space.gradientMatrix';
            E = (1-runData(i).beta)*runData(i).Space.left*runData(i).Space.left' - (1+runData(i).beta)*runData(i).Space.right*runData(i).Space.right';
            Eneg = (1+runData(i).beta)*runData(i).Space.left*runData(i).Space.right';
            Epos = (-1+runData(i).beta)*runData(i).Space.right*runData(i).Space.left';
        end
        for n = 1:numel(k) % loop over wavemodes
            % Residual operator:
            R = runData(i).Space.massMatrix\(2*C + E + exp(-1i*k(n))*Eneg + exp(1i*k(n))*Epos);
            % Modified wavenumbers:
            kMod(:,n) = 1i*eigs(R,J);
            % Make eigenmodes consistent:
            if n > 1
                [~,ids] = min(abs(kMod(:,n) - kMod(:,n-1).'),[],2); % mind the dot! (transpose vs. ctranspose)
                kMod(ids,n) = kMod(:,n);
            end
        end
        % Sort eigenmodes:
        [~,ids] = sort(abs(kMod(:,k==0)));
        kMod = kMod(ids,:);
        % Postprocess & export:
        tic, CFL = runData(i).Time.optimizeCFL(runData(i).Space);
        fprintf(fileID,'# %s, CFL_max = %.12g (%.6g s)\n',class(runData(i).Time),CFL,toc);
        fprintf(fileID,'# %s; J = %d, %d patches, beta = %f\n',runData(i).Space.getName,runData(i).Space.basisCount,runData(i).K,runData(i).beta);
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','j','n','k_exact','k_real','k_imag','z_real','z_imag');
        if all(isinf(runData(i).eigenmodes))
            runData(i).eigenmodes = 1:J; % store all eigenmodes
        end
        runData(i).eigenmodes(runData(i).eigenmodes > J) = []; % throw away inexistent modes
        for j = runData(i).eigenmodes % loop over eigenmodes
            fprintf(fileID,'%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n',[j*ones(1,2*N+1); -N:N; k/J; real(kMod(j,:))/J; imag(kMod(j,:))/J; CFL*imag(kMod(j,:)); CFL*real(kMod(j,:))]);
            fprintf(fileID,'\n');
        end
        fprintf('Run %d of %d completed by worker %d in %.3g s.\n',i,I,get(getCurrentTask(),'ID'),toc)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isnan(runData(i).fig)
            continue
        end
        figure(runData(i).fig)
        hold on
        subplot(2,2,1)
        plot(k/J,real(kMod(runData(i).eigenmodes,:))/J)
        subplot(2,2,3)
        plot(k/J,imag(kMod(runData(i).eigenmodes,:))/J)
        subplot(2,2,[2 4])
        plot(-1i*kMod(runData(i).eigenmodes,:).'*CFL,'.-')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch me
       warning('Run %d of %d did not succeed.\n "%s"',i,I,getReport(me))
    end
    fclose(fileID);
end