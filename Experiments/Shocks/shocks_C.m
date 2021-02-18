clc
clear
%close all

% This script launches all runs of "batch C", which compares DGIGA without
% limiter, with Krivodonova's limiter, and with AFC's limiter.
% No sensors.
% Courant number is fixed.
% Optimal DGIGA bases, of J = 3 or J = 6.

%% Test matrix
dataCell = {
%   run     K       relCFL  basis               	limiter                                                                             									test 
    1       100     .5      DGIGA(1,2)          	Limiter                                                                             									@jiangShu
    2       100     .5      DGIGA(1,2)          	Krivodonova('Stats',true)                                                           									@jiangShu
    3       100     .5      DGIGA_AFC(1,2)      	AFC_2010('Control',1,'Failsafe',3,'Stats',true)                                                             			@jiangShu
	4       100     .5      DGIGA_AFC(1,2)          AFC_2010('Sensor',KXRCF,'Control',1,'Failsafe',3,'Stats',true)                                              			@jiangShu
	5       100     .5      DGIGA_AFC(1,2)          AFC_2010('Sensor',APTVD,'Control',1,'Failsafe',3,'Stats',true)                                              			@jiangShu
										
    6      	100     .5      DGIGA(2,3,1)    		Limiter                                                                             									@jiangShu
    7      	100     .5      DGIGA(2,3,1)    		Krivodonova('Stats',true)                                                           									@jiangShu
    8      	100     .5      DGIGA_AFC(2,3,1)		AFC_2010('Control',1,'Failsafe',3,'Stats',true)                                                             			@jiangShu
    9      	100     .5      DGIGA_AFC(2,3,1)		AFC_2010('Sensor',KXRCF,'Control',1,'Failsafe',3,'Stats',true)                                                      	@jiangShu
    10     	100     .5      DGIGA_AFC(2,3,1)		AFC_2010('Sensor',APTVD,'Control',1,'Failsafe',3,'Stats',true)                                                      	@jiangShu
			
    11	    100     .5      DGIGA(1,2)    			[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@toro1
    12      100     .5      DGIGA(1,2)    			[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@toro1
    13      100     .5      DGIGA_AFC(1,2)			[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@toro1
    14     	100     .5      DGIGA_AFC(1,2)  		[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
    15     	100     .5      DGIGA_AFC(1,2)  		[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro1
		
    16      100     .5      DGIGA(2,3,1)    		[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@toro1
    17      100     .5      DGIGA(2,3,1)    		[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@toro1
    18      100     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@toro1
    19      100     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]	@toro1
    20     	100     .5      DGIGA_AFC(2,3,1)   		[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]	@toro1
				
    21      100     .5      DGIGA(1,2)    			[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@toro2
    22      100     .5      DGIGA(1,2)    			[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@toro2
    23      100     .5      DGIGA_AFC(1,2)			[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@toro2
    24     	100     .5      DGIGA_AFC(1,2)		    [AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    25      100     .5      DGIGA_AFC(1,2)			[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
				
    26      100     .5      DGIGA(2,3,1)    		[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@toro2
    27      100     .5      DGIGA(2,3,1)    		[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@toro2
    28      100     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@toro2
    29      100     .5      DGIGA_AFC(2,3,1)    	[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
    30      100     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @toro2
				
    31      200     .5      DGIGA(1,2)    			[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@woodwardColella
    32      200     .5      DGIGA(1,2)    			[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@woodwardColella
    33      200     .5      DGIGA_AFC(1,2)			[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@woodwardColella
    34     	200     .5      DGIGA_AFC(1,2)   		[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
    35      200     .5      DGIGA_AFC(1,2)			[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella

    36      200     .5      DGIGA(2,3,1)    		[Limiter  EulerP1('Stats',true)   EulerP0('Stats',true)]                   												@woodwardColella
    37      200     .5      DGIGA(2,3,1)    		[Krivodonova('Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   											@woodwardColella
    38      200     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   				@woodwardColella
    39	    200     .5      DGIGA_AFC(2,3,1)    	[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
    40      200     .5      DGIGA_AFC(2,3,1)		[AFC_2010('Sensor',APTVD,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1('Stats',true)   EulerP0('Stats',true)]   @woodwardColella
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
        dataStruct = dataStruct.test(dataStruct,sprintf('%s_%d',batchName,i));
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