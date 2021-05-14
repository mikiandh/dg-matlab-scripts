clc
clear

% This script launches all runs of "batch B", which compares sensors in
% DGSEM, and with Krivodonova's limiter (which seems to be the best of the
% lot), against fine-tuned TVB and each other. Courant number is fixed.
% Degree is 2 or 5.

%% Test matrix
dataCell = {
%   run     K       relCFL  basis       limiter                                                                                     test 
    1       100     .5      DGSEM(2)    Krivodonova('Stats',true)                                                                   @jiangShu
    2       100     .5      DGSEM(2)    Krivodonova('Sensor',KXRCF,'Stats',true)                                                    @jiangShu
    3       100     .5      DGSEM(2)    Krivodonova('Sensor',APTVD,'Stats',true)                                                    @jiangShu
    
    4       100     .5      DGSEM(5)    Krivodonova('Stats',true)                                                                   @jiangShu
    5       100     .5      DGSEM(5)    Krivodonova('Sensor',KXRCF,'Stats',true)                                                    @jiangShu
    6       100     .5      DGSEM(5)    Krivodonova('Sensor',APTVD,'Stats',true)                                                    @jiangShu
    
    7       100     .5      DGSEM(2)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @toro1
    8       100     .5      DGSEM(2)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
    9       100     .5      DGSEM(2)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
    
    10      100     .5      DGSEM(5)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @toro1
    11      100     .5      DGSEM(5)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
    12      100     .5      DGSEM(5)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
    
    13      100     .5      DGSEM(2)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @toro2
    14      100     .5      DGSEM(2)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    15      100     .5      DGSEM(2)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    
    16      100     .5      DGSEM(5)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @toro2
    17      100     .5      DGSEM(5)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    18      100     .5      DGSEM(5)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    
    19      200     .5      DGSEM(2)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @woodwardColella
    20      200     .5      DGSEM(2)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
    21      200     .5      DGSEM(2)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
    
    22      200     .5      DGSEM(5)    [Krivodonova('Stats',true) EulerP1('Stats',true)   EulerP0('Stats',true)]                   @woodwardColella
    23      200     .5      DGSEM(5)    [Krivodonova('Sensor',KXRCF,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
    24      200     .5      DGSEM(5)    [Krivodonova('Sensor',APTVD,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
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