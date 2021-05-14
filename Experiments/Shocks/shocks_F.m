clc
clear

% This script launches all runs of "batch F", which tries to find an
% optimum balance between IGA and DGIGA in the shock-turbulence problem.
% Time-step size is fixed.

%% Test matrix
dataCell = {
%   run     K       relCFL  basis                           limiter                                                                     test 
    1       1       nan     DGIGA_AFC(1198,2,1)             AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true)		@shuOsher
    2       10      nan     DGIGA_AFC(118,2,1)              AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true)		@shuOsher
    3       100     nan     DGIGA_AFC(10,2,1)               AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true)		@shuOsher
    4       200     nan     DGIGA_AFC(4,2,1)                AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true)		@shuOsher
    5       400     nan     DGIGA_AFC(1,2,1)                AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true)		@shuOsher
					
    6       10      nan     DGIGA_AFC(118,2,1)              AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    7       100     nan     DGIGA_AFC(10,2,1)               AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    8       200     nan     DGIGA_AFC(4,2,1)                AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    9       400     nan     DGIGA_AFC(1,2,1)                AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
		
    10      10      nan     DGIGA_AFC(118,2,1)              AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    11      100     nan     DGIGA_AFC(10,2,1)               AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    12      200     nan     DGIGA_AFC(4,2,1)                AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
    13      400     nan     DGIGA_AFC(1,2,1)                AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) 		@shuOsher
};
batchName = mfilename;
ptSkip = repmat(pow2(1:5),1,3); % make the exported solutions lighter by skipping every few points, run-wise

%% Batch setup
[temp{1:size(dataCell,1),1:13}] = deal(nan);
dataStructArray = cell2struct([dataCell temp],...
    {'run','K','relCFL','basis','limiter','test','J',...
    'densityErrorL1','momentumErrorL1','energyErrorL1',...
    'densityTV','momentumTV','energyTV',...
    'densityExactTV','momentumExactTV','energyExactTV',...
    'sensorRatio','limiterRatio','wallClockTime'},2);
I = size(dataStructArray,1);

%% Parallel write setup
c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);
spmd
    cpuFile = (fopen(c.Value));
    fprintf(c.Value,'%s\n',strjoin(fields(dataStructArray),'\t'));
end
    
%% Parallel batch run
parfor i = 1:I
    % Distribute data:
    dataStruct = dataStructArray(i);
    % Try to solve current run:
    try
        dataStruct.J = dataStruct.basis.basisCount;
        dataStruct = dataStruct.test(dataStruct,sprintf('%s_%d',batchName,dataStruct.run),ptSkip(i));
    catch ME
        fprintf('Run %d of %d failed; it will be skipped (worker %d).\n',i,I,get(getCurrentTask(),'ID'))
        continue % do not write it
    end
    % Export to temporary file:
    dataFormats = {
        '%d','%d','%g','%s','%s','%s','%d',...
        '%.12f','%.12f','%.12f',...
        '%.12f','%.12f','%.12f',...
        '%.12f','%.12f','%.12f',...
        '%g','%g','%g'
    };
    fprintf(c.Value,[strjoin(dataFormats,'\t') '\n'],...
        dataStruct.run,...
        dataStruct.K,...
        dataStruct.relCFL,...
        dataStruct.limiter.getInfo,...
        dataStruct.basis.getName,...
        func2str(dataStruct.test),...
        dataStruct.J,...
        dataStruct.densityErrorL1,...
        dataStruct.momentumErrorL1,...
        dataStruct.energyErrorL1,...
        dataStruct.densityTV,...
        dataStruct.momentumTV,...
        dataStruct.energyTV,...
        dataStruct.densityExactTV,...
        dataStruct.momentumExactTV,...
        dataStruct.energyExactTV,...
        dataStruct.sensorRatio,...
        dataStruct.limiterRatio,...
        dataStruct.wallClockTime);
    fprintf('Run %d of %d completed (worker %d).\n',i,I,get(getCurrentTask(),'ID'))
end
clear c % free memory and close temporary files

%% Export to file
try
    % Read temporary files:
    spmd
        cpuTable = readtable(cpuFile,...
            'Delimiter','\t','ReadVariableNames',true); % 'TextType','string' ?
    end
    outputTable = sortrows(vertcat(cpuTable{:}));
    clear cpuTable % free memory
    % Write to output file:
    try
        writetable(outputTable,[batchName '_table.dat'],'Delimiter','\t') % 'QuoteStrings', true ?
        spmd
            delete(cpuFile) % delete them
        end
    catch me
        warning('Final export to file failed.')
    end
catch me
    warning('Could not read temporary files.')
end