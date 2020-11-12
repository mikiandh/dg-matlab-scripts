% This script computes and writes into files the estimated cost (in FLOPs)
% of solving advection using RK3 and some FR scheme.
%
clear, clc, close all

%% File name root
fileRoot = 'flops';

%% Methods
bases = [
    % optimal FR                       % DGSEM via FR   % DGSEM
    FR({'eta',0.0227795117943271},2)   FR({'eta',0},2)  DGSEM(2)
    FR({'eta',0.0351005869100722},3)   FR({'eta',0},3)  DGSEM(3)
    FR({'eta',0.047091694344944},4)    FR({'eta',0},4)  DGSEM(4)
    FR({'eta',0.0580283806157087},5)   FR({'eta',0},5)  DGSEM(5)
    FR({'eta',0.0786884021016807},7)   FR({'eta',0},7)  DGSEM(7)
    FR({'eta',0.102887003589808},10)   FR({'eta',0},10) DGSEM(10)
    FR({'eta',0.128444040756296},14)   FR({'eta',0},14) DGSEM(14)
    FR({'eta',0.152203674556189},19)   FR({'eta',0},19) DGSEM(19)
    ]';

%% Test conditions
waveRange0 = logspace(1,5,25)'; % range of wavemodes to resolve (seed)
timeRange = 1; % simulated time span, in domain lengths crossed
solver = SSP_RK3;

%% Parallel loop over basis pairs:
parfor j = 1:size(bases,2)
    tbl = table;
    for i = 1:size(bases,1)
        % Get basis properties:
        J = bases(i,j).basisCount;
        kf = bases(i,j).getResolvingWavenumber;
        cfl = solver.optimizeCFL(bases(i,j)); %#ok<PFBNS>
        % Deduce (integer) K to achieve sufficient resolution:
        K = 2*pi*waveRange0/kf;
        K = ceil(K);
        waveRange = K*kf/(2*pi); % actual well-resolved wavemode ratios
        % Deduce integer number of time-steps:
        N = timeRange*K/solver.optimizeCFL(bases(i,j));
        N = ceil(N); % last step is smaller but costs the same
        if i == 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flops = 2*K*J^2 + 23*J*K + 8*J + 6*K + 5; %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flops = 2*K*J^2 + 25*J*K + 8*J + 8*K + 5; %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        flops = flops.*N;
        if i == 1
            tbl{:,1:4} = [waveRange K N flops];
        else
            tbl{:,3+i} = (tbl{:,4} - flops)./flops;
            hold on
        end
        plot(waveRange,flops,'x-','DisplayName',bases(i,j).getName)
        hold off
    end
    tbl.Properties.VariableNames = {'scaleRatio','K','steps','flops','relFlopsEta0','relFlopsDGSEM'};
    % Print to file:
    fileName = sprintf('%s_%s',fileRoot,strjoin(regexp(bases(1,j).getName,'\d+|^[A-Z]+|\((\w+)\)','match'),'_'));
    writetable(tbl,[fileName '.dat'],'Delimiter','\t');
%     % Export figure:
%     legend('Location','Best')
%     saveas(gcf,fileName,'png')
end