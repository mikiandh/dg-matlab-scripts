clc
clear
%close all

% This script launches all runs of "batch A", which compares all inter-cell
% limiters (all vs. all) in DGSEM. No sensor is used. Courant number is
% fixed. Degree is 2 or 5.

%% Batch parameters
K = { % number of patches/elements
    100
    };
relCFL = { % Courant number, relative to max. linearly stable
    .5
    };
bases = {
    DGSEM(2)
    DGSEM(5)
    };
limiters = {
    [Limiter EulerP1 EulerP0]
    [TVB EulerP1 EulerP0]
    [BDF EulerP1 EulerP0]
    [BSB EulerP1 EulerP0]
    [Krivodonova EulerP1 EulerP0]
    [WENO EulerP1 EulerP0]
};
tests = {
    @jiangShu
    @toro1
    @toro2
    @woodwardColella
};
batchName = mfilename;

%% Batch setup
dataTable = cell2table(...
    assembleTestTuple(limiters,bases,tests,relCFL,K),...
    'VariableNames',{'run','limiter','basis','test','relCFL','K'});
dataTable.J = arrayfun(@(x) x.basisCount,dataTable.basis);
I = size(dataTable,1);
disp(dataTable)

% Preallocate extra columns:
dataTable.densityErrorL1(:) = nan;
dataTable.momentumErrorL1(:) = nan;
dataTable.energyErrorL1(:) = nan;
dataTable.densityTV(:) = nan;
dataTable.momentumTV(:) = nan;
dataTable.energyTV(:) = nan;
dataTable.densityExactTV(:) = nan;
dataTable.momentumExactTV(:) = nan;
dataTable.energyExactTV(:) = nan;
dataTable.troubledDofs(:) = nan;
dataTable.limitedDofs(:) = nan;
dataTable.wallClockTime(:) = nan;

%% Parallel write setup
c = parallel.pool.Constant(@() fopen(tempname(pwd),'wt'),@fclose);
spmd
    cpuFile = (fopen(c.Value));
    fprintf(c.Value,'%s\n',...
        strjoin(dataTable.Properties.VariableNames,' \t'));
end
    
%% Parallel batch run
parfor i = 1:I
    % Distribute data:
    dataRow = dataTable(i,:);
    % Try to solve current run:
    try
        test = dataRow.test{1};
        dataRow = test(dataRow,sprintf('%s_%d',batchName,i));
    catch
        fprintf('Run %d of %d failed; it will be skipped (worker %d).\n',i,I,get(getCurrentTask(),'ID'))
        continue % do not write it
    end
    % Export to temporary file:
    limiters = dataRow{1,2};
    fprintf(c.Value,...
        '%d\t%s\t%s\t%s\t%.3f\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%d\t%d\t%.3f\n',...
        dataRow{1,1},limiters.getInfo,dataRow{1,3}.getName,func2str(dataRow{1,4}{1}),dataRow{1,5:end});
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