clc
clear
%close all

% This script launches all runs of 'batch G'.
% It tries to answer the question: which limiters (if any) preserve order
% of accuracy? Do they, really (in practice)? 
% Courant number fixed.

%% Test matrix
dataCell = {
%   run     K       relCFL      basis               limiter                                                                 test
    1       24      .5          DGIGA(1,1)          TVB('Stats',true)                                                       @gaussHump
    2       36      .5          DGIGA(1,1)          TVB('Stats',true)                                                       @gaussHump
    3       48      .5          DGIGA(1,1)          TVB('Stats',true)                                                       @gaussHump

    4       16      .5          DGIGA(1,2)          TVB('Stats',true)                                                       @gaussHump
    5       24      .5          DGIGA(1,2)          TVB('Stats',true)                                                       @gaussHump
    6       32      .5          DGIGA(1,2)          TVB('Stats',true)                                                       @gaussHump
    
    7       12      .5          DGIGA(1,3)          TVB('Stats',true)                                                       @gaussHump
    8       18      .5          DGIGA(1,3)          TVB('Stats',true)                                                       @gaussHump
    9       24      .5          DGIGA(1,3)          TVB('Stats',true)                                                       @gaussHump
    
    10      24      .5          DGIGA(1,1)          Krivodonova('Stats',true)                                               @gaussHump
    11      36      .5          DGIGA(1,1)          Krivodonova('Stats',true)                                               @gaussHump
    12      48      .5          DGIGA(1,1)          Krivodonova('Stats',true)                                               @gaussHump
    
    13      16      .5          DGIGA(1,2)          Krivodonova('Stats',true)                                               @gaussHump
    14      24      .5          DGIGA(1,2)          Krivodonova('Stats',true)                                               @gaussHump
    15      32      .5          DGIGA(1,2)          Krivodonova('Stats',true)                                               @gaussHump
    
    16      12      .5          DGIGA(1,3)          Krivodonova('Stats',true)                                               @gaussHump
    17      18      .5          DGIGA(1,3)          Krivodonova('Stats',true)                                               @gaussHump
    18      24      .5          DGIGA(1,3)          Krivodonova('Stats',true)                                               @gaussHump
    
    19      24      .5          DGIGA_AFC(1,1)      AFC_2010('Stats',true)                                                  @gaussHump
    20      36      .5          DGIGA_AFC(1,1)      AFC_2010('Stats',true)                                                  @gaussHump
    21      48      .5          DGIGA_AFC(1,1)      AFC_2010('Stats',true)                                                  @gaussHump
    
    22      16      .5          DGIGA_AFC(1,2)      AFC_2010('Stats',true)                                                  @gaussHump
    23      24      .5          DGIGA_AFC(1,2)      AFC_2010('Stats',true)                                                  @gaussHump
    24      32      .5          DGIGA_AFC(1,2)      AFC_2010('Stats',true)                                                  @gaussHump
    
    25      12      .5          DGIGA_AFC(1,3)      AFC_2010('Stats',true)                                                  @gaussHump
    26      18      .5          DGIGA_AFC(1,3)      AFC_2010('Stats',true)                                                  @gaussHump
    27      24      .5          DGIGA_AFC(1,3)      AFC_2010('Stats',true)                                                  @gaussHump

    28      24      .5          DGSEM(1)            Krivodonova('Stats',true)                                               @gaussHump
    29      36      .5          DGSEM(1)            Krivodonova('Stats',true)                                               @gaussHump
    30      48      .5          DGSEM(1)            Krivodonova('Stats',true)                                               @gaussHump
    
    31      16      .5          DGSEM(2)            Krivodonova('Stats',true)                                               @gaussHump
    32      24      .5          DGSEM(2)            Krivodonova('Stats',true)                                               @gaussHump
    33      32      .5          DGSEM(2)            Krivodonova('Stats',true)                                               @gaussHump
    
    34      12      .5          DGSEM(3)            Krivodonova('Stats',true)                                               @gaussHump
    35      18      .5          DGSEM(3)            Krivodonova('Stats',true)                                               @gaussHump
    36      24      .5          DGSEM(3)            Krivodonova('Stats',true)                                               @gaussHump    

    37      1       .5          DGIGA_AFC(47,1)     AFC_2010('Stats',true)                                                  @gaussHump
    38      1       .5          DGIGA_AFC(71,1)     AFC_2010('Stats',true)                                                  @gaussHump
    39      1       .5          DGIGA_AFC(95,1)     AFC_2010('Stats',true)                                                  @gaussHump
    
    40      1       .5          DGIGA_AFC(46,2)     AFC_2010('Stats',true)                                                  @gaussHump
    41      1       .5          DGIGA_AFC(70,2)     AFC_2010('Stats',true)                                                  @gaussHump
    42      1       .5          DGIGA_AFC(94,2)     AFC_2010('Stats',true)                                                  @gaussHump
    
    43      1       .5          DGIGA_AFC(45,3)     AFC_2010('Stats',true)                                                  @gaussHump
    44      1       .5          DGIGA_AFC(69,3)     AFC_2010('Stats',true)                                                  @gaussHump
    45      1       .5          DGIGA_AFC(93,3)     AFC_2010('Stats',true)                                                  @gaussHump    

    46      300     .5          DGIGA(1,1)          TVB('Sensor',KXRCF,'Stats',true)                                        @shuOsher
    47      300     .5          DGIGA(1,2)          TVB('Sensor',KXRCF,'Stats',true)                                        @shuOsher
    48      300     .5          DGIGA(1,3)          TVB('Sensor',KXRCF,'Stats',true)                                        @shuOsher
    
    49      300     .5          DGIGA(1,1)          Krivodonova('Sensor',KXRCF,'Stats',true)                                @shuOsher
    50      300     .5          DGIGA(1,2)          Krivodonova('Sensor',KXRCF,'Stats',true)                                @shuOsher
    51      300     .5          DGIGA(1,3)          Krivodonova('Sensor',KXRCF,'Stats',true)                                @shuOsher
    
    52      300     .5          DGIGA_AFC(1,1)      AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)    @shuOsher
    53      300     .5          DGIGA_AFC(1,2)      AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)    @shuOsher
    54      300     .5          DGIGA_AFC(1,3)      AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)    @shuOsher    
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

%% Compute convergence orders
outputTable.order(:) = nan;
for i = [1:3:45; 2:3:45; 3:3:45]
    j = find(any(outputTable.run' == i));
    logError = log(outputTable(j,:).densityErrorL1);
    logElems = log(outputTable(j,:).K);
    outputTable.order(j(2:end)) = diff(logError)./diff(logElems);
end