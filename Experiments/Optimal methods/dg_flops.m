% This script computes and writes into files the estimated cost (in FLOPs)
% of solving advection using RK3 and some DG scheme.
%
clear, clc, close all

%% File name root
fileRoot = 'flops';

%% Methods
bases = [
    DGSEM(0) DGSEM(1) DGSEM(2) DGSEM(3) DGSEM(4) DGSEM(5) DGSEM(7) DGSEM(10) DGSEM(14) DGSEM(19)     % baselines
    DGSEM(0) DGSEM(0) DGSEM(0) DGSEM(0) DGSEM(0) DGSEM(0) DGSEM(0) DGSEM(0)  DGSEM(0)  DGSEM(0)      % baselines' baseline
    ];

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        flops = 2*K.*J^2 + 23*J*K + 8*J + 6*K + 5;                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        flops = flops.*N;
        if i == 1
            tbl{:,1:4} = [waveRange K N flops];
        else
            tbl{:,5} = (tbl{:,4} - flops)./flops;
            hold on
        end
        plot(waveRange,flops,'x-','DisplayName',bases(i,j).getName)
        hold off
    end
    tbl.Properties.VariableNames = {'scaleRatio','K','steps','flops','relFlops'};
    % Print to file:
    fileName = sprintf('%s_%s',fileRoot,strjoin(regexp(bases(1,j).getName,'\d+|^[A-Z]+|\((\w+)\)','match'),'_'));
    writetable(tbl,[fileName '.dat'],'Delimiter','\t');
%     % Export figure:
%     legend('Location','Best')
%     saveas(gcf,fileName,'png')
end