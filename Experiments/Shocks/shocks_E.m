clc
clear
%close all

% This script launches all runs of "batch E", which compares the best
% limiting solutions for the three DG variants, in the shock-turbulence
% problem.
% Courant number is fixed.

%% Test matrix
dataCell = {
%   run     K       relCFL  basis                                   limiter                                                                             											test 
    1       400     .5      DGSEM(2)                                [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    2       400     .5      FR({'eta',0.0227795117943271},2)		[Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    3       400     .5      DGIGA(1,2,1)                            [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
	4       400     .5      DGIGA_AFC(1,2,1)                        [AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]  	@shuOsher
	
    5       400     .5      DGSEM(2)                                [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    6       400     .5      FR({'eta',0.0227795117943271},2)		[Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    7       400     .5      DGIGA(1,2,1)                            [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
	8       400     .5      DGIGA_AFC(1,2,1)                        [AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    9       400     .5      DGSEM(2)                                [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    10      400     .5      FR({'eta',0.0227795117943271},2)        [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    11      400     .5      DGIGA(1,2,1)                            [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
	12      400     .5      DGIGA_AFC(1,2,1)                        [AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    13      200     .5      DGSEM(5)                                [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    14      200     .5      FR({'eta',0.0580283806157087},5)		[Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    15      200     .5      DGIGA(2,3,1)                            [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
	16      200     .5      DGIGA_AFC(2,3,1)                        [AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    17      200     .5      DGSEM(5)                                [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    18      200     .5      FR({'eta',0.0580283806157087},5)		[Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    19      200     .5      DGIGA(2,3,1)                            [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    20      200     .5      DGIGA_AFC(2,3,1)                        [AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    21      200     .5      DGSEM(5)                                [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    22      200     .5      FR({'eta',0.0580283806157087},5)        [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    23      200     .5      DGIGA(2,3,1)                            [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    24      200     .5      DGIGA_AFC(2,3,1)                        [AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
									
    25      60      .5      DGSEM(19)                               [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    26      60      .5      FR({'eta',0.152203674556189},19) 		[Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    27      60      .5      DGIGA(2,14,9)                           [Krivodonova('Sensor',Sensor,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]         								@shuOsher
    28      60      .5      DGIGA_AFC(2,14,9)                       [AFC_2010('Sensor',Sensor,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    29      60      .5      DGSEM(19)                               [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    30      60      .5      FR({'eta',0.152203674556189},19) 		[Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    31      60      .5      DGIGA(2,14,9)                           [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    32      60      .5      DGIGA_AFC(2,14,9)                       [AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
								
    33      60      .5      DGSEM(19)                               [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    34      60      .5      FR({'eta',0.152203674556189},19)        [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    35      60      .5      DGIGA(2,14,9)                           [Krivodonova('Sensor',APTVD,'Stats',true) EulerP1('Stats',true) EulerP0('Stats',true)]          								@shuOsher
    36      60      .5      DGIGA_AFC(2,14,9)                       [AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true) EulerP1_step('Stats',true) EulerP0_step('Stats',true)]	@shuOsher
};
batchName = mfilename;

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
        dataStruct = dataStruct.test(dataStruct,sprintf('%s_%d',batchName,dataStruct.run));
    catch
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