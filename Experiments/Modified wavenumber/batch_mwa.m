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
% DGSEM, batch 5
    'Filename'              'Time'     'Space'         'K'      'beta'
    'mwa_dgsem_1.dat'       SSP_RK3    DGIGA(1,1,1)    60       1
};

%% Setup
tblIn = cell2table(inputData(2:end,2:end),...
    'VariableNames',inputData(1,2:end),...
    'RowNames',inputData(2:end,1));
runData = table2struct(tblIn);
if any(any(~cellfun(@isscalar,inputData(2:end,2:end)))) % safety check
    error('Please, specify one run per file only.')
end
tblIn.Time = arrayfun(@(x) class(x),tblIn.Time,'UniformOutput',false);
tblIn.Space = arrayfun(@(x) x.getName,tblIn.Space,'UniformOutput',false);
disp(tblIn) % pretty print
I = numel(runData);

%% Parallel run
for i = 1:I
    fileID = fopen(tblIn(i,:).Row{1},'w');
    try
        % Wavenumbers and stuff:
        J = runData(i).Space.basisCount; % number of eigenmodes (i.e. basis functions per patch)
        N = round(runData(i).K*J/2); % number of wavemodes
        k = 2*pi/runData(i).K*(1:N); % all positive (nondimensional) wavenumbers
        kMod = nan(J,N); % corresponding modified wavenumbers
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
        for n = 1:N % loop over wavemodes
            % Residual operator:
            R = runData(i).Space.massMatrix\(2*C + E + exp(-1i*k(n))*Eneg + exp(1i*k(n))*Epos);
            % Modified wavenumbers:
            kMod(:,n) = 1i*eigs(R,J,'smallestabs'); % places physical mode first
            % Sort them:
            if n > 1
                [~,ids] = min(abs(kMod(:,n) - kMod(:,n-1).'),[],2); % mind the dot! (transpose vs. ctranspose)
                kMod(ids,n) = kMod(:,n); % modes [1 2 3] assigned to <ids>
            end
        end
        %%%
        figure
        subplot(2,2,1)
        plot(k/J,real(kMod(2:end,:))/J,'--',k/J,real(kMod(1,:))/J,'-')
        subplot(2,2,2)
        plot(k/J,imag(kMod(2:end,:))/J,'--',k/J,imag(kMod(1,:))/J,'-')
        subplot(2,2,[3 4])
        plot(-1i*kMod.','-')
        %%%
        % Postprocess & export:
        tic, CFL = runData(i).Time.optimizeCFL(runData(i).Space);
        fprintf(fileID,'# %s, CFL_max = %.12g (%.6g s)\n',class(runData(i).Time),CFL,toc);
        fprintf(fileID,'# %s; J = %d, %d patches, beta = %f\n',runData(i).Space.getName,runData(i).Space.basisCount,runData(i).K,runData(i).beta);
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n','j','n','k_real','k_imag','k_exact');
        for j = 1:J % loop over eigenmodes
            fprintf(fileID,'%d\t%d\t%.12g\t%.12g\t%.12g\n',[j*ones(1,N); 1:N; real(kMod(j,:))/J; imag(kMod(j,:))/J; k/J]);
            fprintf(fileID,'\n'); % blank line
        end
        fprintf('Run %d of %d completed by worker %d in %.3g s.\n',i,I,get(getCurrentTask(),'ID'),toc)
    catch me
       warning('Run %d of %d did not succeed.\n "%s"',i,I,getReport(me))
    end
    fclose(fileID);
end