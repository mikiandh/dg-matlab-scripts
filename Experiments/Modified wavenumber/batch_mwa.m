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
    'Filename'                  'Time'     'Space'         	'K'      'beta'     'eigenmodes'     'fig'
% % DGSEM
%     'mwa_fv_upwind.dat'         SSP_RK3    DGSEM(0)        	60       1          inf              1
%     'mwa_fv_center.dat'         SSP_RK3    DGSEM(0)        	60       0          inf              1
% %-------------------------------------------------------------------------------------------------------
%     'mwa_dgsem_2_all.dat'       SSP_RK3    DGSEM(2)        	60       1          inf              2
%     'mwa_dgsem_3_all.dat'       SSP_RK3    DGSEM(3)        	60       1          inf              3
%     'mwa_dgsem_4_all.dat'       SSP_RK3    DGSEM(4)        	60       1          inf              4
% %-------------------------------------------------------------------------------------------------------
%     'mwa_dgsem_0.dat'           SSP_RK3    DGSEM(0)        	60       1          1                5
%     'mwa_dgsem_1.dat'           SSP_RK3    DGSEM(1)        	60       1          1                5
%     'mwa_dgsem_2.dat'           SSP_RK3    DGSEM(2)        	60       1          1                5
%     'mwa_dgsem_3.dat'           SSP_RK3    DGSEM(3)        	60       1          1                5
%     'mwa_dgsem_4.dat'           SSP_RK3    DGSEM(4)        	60       1          1                5
%     'mwa_dgsem_5.dat'           SSP_RK3    DGSEM(5)        	60       1          1                5 
% %-------------------------------------------------------------------------------------------------------
%     'mwa_dgsem_9.dat'           SSP_RK3    DGSEM(9)        	60       1          1                6
%     'mwa_dgsem_18.dat'          SSP_RK3    DGSEM(18)       	60       1          1                6
%     'mwa_dgsem_33.dat'          SSP_RK3    DGSEM(33)       	60       1          1                6
%     'mwa_dgsem_63.dat'          SSP_RK3    DGSEM(63)       	60       1          1                6
%     'mwa_dgsem_119.dat'         SSP_RK3    DGSEM(119)      	60       1          1                6
% % FR	
%     'mwa_fr_2_min_all.dat'      SSP_RK3    FR('min',2)     	60       1          inf              7
%     'mwa_fr_2_dg_all.dat'       SSP_RK3    FR('DG',2)     	60       1          inf              7
%     'mwa_fr_2_ga_all.dat'       SSP_RK3    FR('Ga',2)      	60       1          inf              7
%     'mwa_fr_2_lulo_all.dat'     SSP_RK3    FR('LumpLo',2)  	60       1          inf              7
%     'mwa_fr_2_max_all.dat'      SSP_RK3    FR('max',2)     	60       1          inf              7
%     'mwa_fr_3_min_all.dat'      SSP_RK3    FR('min',3)     	60       1          inf              8
%     'mwa_fr_3_dg_all.dat'       SSP_RK3    FR('DG',3)     	60       1          inf              8
%     'mwa_fr_3_ga_all.dat'       SSP_RK3    FR('Ga',3)      	60       1          inf              8
%     'mwa_fr_3_lulo_all.dat'     SSP_RK3    FR('LumpLo',3)  	60       1          inf              8
%     'mwa_fr_3_max_all.dat'      SSP_RK3    FR('max',3)     	60       1          inf              8
% %-------------------------------------------------------------------------------------------------------
%     'mwa_fr_2_min.dat'          SSP_RK3    FR('min',2)     	60       1          1                9
%     'mwa_fr_3_min.dat'          SSP_RK3    FR('min',3)     	60       1          1                9
%     'mwa_fr_4_min.dat'          SSP_RK3    FR('min',4)     	60       1          1                9
%     'mwa_fr_5_min.dat'          SSP_RK3    FR('min',5)     	60       1          1                9
%     'mwa_fr_8_min.dat'          SSP_RK3    FR('min',8)     	60       1          1                9
%     'mwa_fr_30_min.dat'         SSP_RK3    FR('min',30)    	60       1          1                9
%     'mwa_fr_119_min.dat'        SSP_RK3    FR('min',119)   	60       1          1                9
% %-------------------------------------------------------------------------------------------------------
%     'mwa_fr_2_dg.dat'           SSP_RK3    FR('DG',2)      	60       1          1                10
%     'mwa_fr_3_dg.dat'           SSP_RK3    FR('DG',3)      	60       1          1                10
%     'mwa_fr_4_dg.dat'           SSP_RK3    FR('DG',4)      	60       1          1                10
%     'mwa_fr_5_dg.dat'           SSP_RK3    FR('DG',5)      	60       1          1                10
%     'mwa_fr_8_dg.dat'           SSP_RK3    FR('DG',8)      	60       1          1                10
%     'mwa_fr_30_dg.dat'          SSP_RK3    FR('DG',30)     	60       1          1                10
%     'mwa_fr_119_dg.dat'         SSP_RK3    FR('DG',119)    	60       1          1                10
% %-------------------------------------------------------------------------------------------------------
%     'mwa_fr_2_ga.dat'           SSP_RK3    FR('Ga',2)      	60       1          1                11
%     'mwa_fr_3_ga.dat'           SSP_RK3    FR('Ga',3)      	60       1          1                11
%     'mwa_fr_4_ga.dat'           SSP_RK3    FR('Ga',4)      	60       1          1                11
%     'mwa_fr_5_ga.dat'           SSP_RK3    FR('Ga',5)      	60       1          1                11
%     'mwa_fr_8_ga.dat'           SSP_RK3    FR('Ga',8)      	60       1          1                11
%     'mwa_fr_30_ga.dat'          SSP_RK3    FR('Ga',30)     	60       1          1                11
%     'mwa_fr_119_ga.dat'         SSP_RK3    FR('Ga',119)    	60       1          1                11
% %-------------------------------------------------------------------------------------------------------
%     'mwa_fr_2_lulo.dat'         SSP_RK3    FR('LumpLo',2)  	60       1          1                12
%     'mwa_fr_3_lulo.dat'         SSP_RK3    FR('LumpLo',3)  	60       1          1                12
%     'mwa_fr_4_lulo.dat'         SSP_RK3    FR('LumpLo',4)  	60       1          1                12
%     'mwa_fr_5_lulo.dat'         SSP_RK3    FR('LumpLo',5)  	60       1          1                12
%     'mwa_fr_8_lulo.dat'         SSP_RK3    FR('LumpLo',8)  	60       1          1                12
%     'mwa_fr_30_lulo.dat'        SSP_RK3    FR('LumpLo',30) 	60       1          1                12
%     'mwa_fr_119_lulo.dat'       SSP_RK3    FR('LumpLo',119)	60       1          1                12
% %-------------------------------------------------------------------------------------------------------
%     'mwa_fr_2_max.dat'          SSP_RK3    FR('max',2)   	60       1          1                13
%     'mwa_fr_3_max.dat'          SSP_RK3    FR('max',3)   	60       1          1                13
%     'mwa_fr_4_max.dat'          SSP_RK3    FR('max',4)   	60       1          1                13
%     'mwa_fr_5_max.dat'          SSP_RK3    FR('max',5)   	60       1          1                13
%     'mwa_fr_8_max.dat'          SSP_RK3    FR('max',8)   	60       1          1                13
%     'mwa_fr_30_max.dat'         SSP_RK3    FR('max',30)  	60       1          1                13
%     'mwa_fr_119_max.dat'        SSP_RK3    FR('max',119) 	60       1          1                13
% DGIGA
    'mwa_dgiga_1_2.dat'         SSP_RK3   DGIGA(1,2)      	60       1          1                17
    'mwa_dgiga_1_3.dat'         SSP_RK3   DGIGA(1,3)        60       1          1                17
    'mwa_dgiga_1_4.dat'         SSP_RK3   DGIGA(1,4)     	60       1          1                17
    'mwa_dgiga_1_5.dat'         SSP_RK3   DGIGA(1,5)     	60       1          1                17
    
    'mwa_dgiga_2_2_0.dat'       SSP_RK3    DGIGA(2,2,0)   	60       1          1                18
    'mwa_dgiga_2_2.dat'         SSP_RK3    DGIGA(2,2)   	60       1          1                19
    'mwa_dgiga_2_3_0.dat'       SSP_RK3    DGIGA(2,3,0)   	60       1          1                18
    'mwa_dgiga_2_3.dat'         SSP_RK3    DGIGA(2,3)   	60       1          1                19
    'mwa_dgiga_2_4_0.dat'       SSP_RK3    DGIGA(2,4,0)   	60       1          1                18
	'mwa_dgiga_2_4_1.dat'       SSP_RK3    DGIGA(2,4,1)   	60       1          1                26
	'mwa_dgiga_2_4_2.dat'       SSP_RK3    DGIGA(2,4,2)   	60       1          1                26
    'mwa_dgiga_2_4.dat'         SSP_RK3    DGIGA(2,4)   	60       1          1                19
    'mwa_dgiga_2_5_0.dat'       SSP_RK3    DGIGA(2,5,0)   	60       1          1                18
    'mwa_dgiga_2_5.dat'         SSP_RK3    DGIGA(2,5)   	60       1          1                19
    
    'mwa_dgiga_4_2_0.dat'       SSP_RK3    DGIGA(4,2,0)   	60       1          1                20
    'mwa_dgiga_4_2.dat'         SSP_RK3    DGIGA(4,2)   	60       1          1                21
    'mwa_dgiga_4_3_0.dat'       SSP_RK3    DGIGA(4,3,0)   	60       1          1                20
    'mwa_dgiga_4_3.dat'         SSP_RK3    DGIGA(4,3)   	60       1          1                21
    'mwa_dgiga_4_4_0.dat'       SSP_RK3    DGIGA(4,4,0)   	60       1          1                20
	'mwa_dgiga_4_4_1.dat'       SSP_RK3    DGIGA(4,4,1)   	60       1          1                27
	'mwa_dgiga_4_4_2.dat'       SSP_RK3    DGIGA(4,4,2)   	60       1          1                27
    'mwa_dgiga_4_4.dat'         SSP_RK3    DGIGA(4,4)   	60       1          1                21
    'mwa_dgiga_4_5_0.dat'       SSP_RK3    DGIGA(4,5,0)   	60       1          1:3              20
    'mwa_dgiga_4_5.dat'         SSP_RK3    DGIGA(4,5)   	60       1          1                21
    
    'mwa_dgiga_8_2_0.dat'       SSP_RK3    DGIGA(8,2,0)   	60       1          1                22
    'mwa_dgiga_8_2.dat'         SSP_RK3    DGIGA(8,2)   	60       1          1                23
    'mwa_dgiga_8_3_0.dat'       SSP_RK3    DGIGA(8,3,0)   	60       1          1:3              22
    'mwa_dgiga_8_3.dat'         SSP_RK3    DGIGA(8,3)   	60       1          1                23
    'mwa_dgiga_8_4_0.dat'       SSP_RK3    DGIGA(8,4,0)   	60       1          1:5              22
	'mwa_dgiga_8_4_1.dat'       SSP_RK3    DGIGA(8,4,1)   	60       1          1                28
	'mwa_dgiga_8_4_2.dat'       SSP_RK3    DGIGA(8,4,2)   	60       1          1                28
    'mwa_dgiga_8_4.dat'         SSP_RK3    DGIGA(8,4)   	60       1          1                23
    'mwa_dgiga_8_5_0.dat'       SSP_RK3    DGIGA(8,5,0)   	60       1          1:7              22
    'mwa_dgiga_8_5.dat'         SSP_RK3    DGIGA(8,5)   	60       1          1                23
	
	'mwa_dgiga_16_2.dat'        SSP_RK3    DGIGA(16,2)   	60       1          1                24
    'mwa_dgiga_16_3.dat'        SSP_RK3    DGIGA(16,3)   	60       1          1                24
    'mwa_dgiga_16_4.dat'        SSP_RK3    DGIGA(16,4)   	60       1          1                24
    'mwa_dgiga_16_5.dat'        SSP_RK3    DGIGA(16,5)   	60       1          1                24
	
	'mwa_dgiga_32_2.dat'        SSP_RK3    DGIGA(32,2)   	60       1          1:3              25
    'mwa_dgiga_32_3.dat'        SSP_RK3    DGIGA(32,3)   	60       1          1                25
    'mwa_dgiga_32_4.dat'        SSP_RK3    DGIGA(32,4)   	60       1          1                25
    'mwa_dgiga_32_5.dat'        SSP_RK3    DGIGA(32,5)   	60       1          1                25
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
        if isempty(getCurrentTask()) && ~isnan(runData(i).fig)
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
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch me
       warning('Run %d of %d did not succeed.\n "%s"',i,I,getReport(me))
    end
    fclose(fileID);
end