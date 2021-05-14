clc
clear

% This script launches all runs of "batch D", which tries to make IGA-AFC
% work.
% Time-step size is fixed.
% Single patch of p = 1,2,3.
% Smoothnesses: C^{0} to C^{p-1}.
% Knot spans: as many as necessary to get ~300 DOFs.

%% Test matrix
dataCell = {
%   run     K       relCFL  basis       				limiter                                                    	test 
    1       1       nan     DGIGA_AFC(299,1,0)  		AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
    2       1       nan     DGIGA_AFC(150,2,0)  		AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
    3       1       nan     DGIGA_AFC(298,2,1)  		AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
    4       1       nan     DGIGA_AFC(100,3,0)  		AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
    5       1       nan     DGIGA_AFC(149,3,1)  		AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
	6       1       nan     DGIGA_nodal_AFC(297,3,2)  	AFC_2010('Control',1,'Failsafe',3,'Stats',true)            	@jiangShu
		
    7       1       nan     DGIGA_AFC(299,1,0)  		AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
    8       1       nan     DGIGA_AFC(150,2,0)  		AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
    9       1       nan     DGIGA_AFC(298,2,1)  		AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
    10      1       nan     DGIGA_AFC(100,3,0)  		AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
    11      1       nan     DGIGA_AFC(149,3,1)  		AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
    12      1       nan     DGIGA_nodal_AFC(297,3,2)  	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro1
			
	13      1       nan     DGIGA_AFC(299,1,0) 			AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
    14      1       nan     DGIGA_AFC(150,2,0) 			AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
    15      1       nan     DGIGA_AFC(298,2,1) 			AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
    16      1       nan     DGIGA_AFC(100,3,0) 			AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
    17      1       nan     DGIGA_AFC(149,3,1) 			AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
	18      1       nan     DGIGA_nodal_AFC(297,3,2) 	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@toro2
		
	19      1       nan     DGIGA_nodal_AFC(599,1,0)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella
    20      1       nan     DGIGA_nodal_AFC(300,2,0)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella
    21      1       nan     DGIGA_nodal_AFC(598,2,1)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella
    22      1       nan     DGIGA_nodal_AFC(200,3,0)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella
    23      1       nan     DGIGA_nodal_AFC(299,3,1)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella
    24      1       nan     DGIGA_nodal_AFC(597,3,2)	AFC_2010('Control',[1 3 2],'Failsafe',3,'Stats',true)      	@woodwardColella	
};
for basis = [dataCell{[dataCell{:,1}] >= 13 & [dataCell{:,1}] <= 18,4}]
    basis.diffusionFun = @DGIGA_AFC.diffusionRobust;
end
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
    fprintf(c.Value,'%s\n',...
        strjoin(fields(dataStructArray),'\t'));
end
    
%% Parallel batch run
parfor i = 1:I
    % Distribute data:
    dataStruct = dataStructArray(i);
    % Try to solve current run:
    try
        dataStruct.J = dataStruct.basis.basisCount;
        dataStruct = dataStruct.test(dataStruct,sprintf('%s_%d',batchName,dataStruct.run),2);
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