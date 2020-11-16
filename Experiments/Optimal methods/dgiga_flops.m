% This script computes and writes into files the estimated cost (in FLOPs)
% of solving advection using RK3 and some DGIGA scheme.
%
clear, clc, close all

%% File name root
fileRoot = 'flops';

%% Methods
bases = [
    DGIGA(1,2,1) DGIGA(1,3,2) DGIGA(1,4,3) DGIGA(2,3,1) DGIGA(2,5,3) DGIGA(2,7,4)  DGIGA(2,10,6)  DGIGA(2,14,9)  % optima
    DGIGA(1,2,1) DGIGA(1,3,2) DGIGA(1,4,3) DGIGA(1,5,4) DGIGA(1,7,6) DGIGA(1,10,9) DGIGA(1,14,13) DGIGA(1,19,18) % Bernstein
    DGSEM(2)     DGSEM(3)     DGSEM(4)     DGSEM(5)     DGSEM(7)     DGSEM(10)     DGSEM(14)      DGSEM(19)      % baselines
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
        if i == 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flops = 6*K.*J^2 + 45*J*K + 24*J + 18*K + 15; %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            k = bases(i,j).nonzeroSpanCount;
            p = bases(i,j).degree;
            s = bases(i,j).smoothness;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flops = 24*J + 36*K - 6*K*p + 12*K*s + 6*J^2*K - 6*K*p^2 + 12*K*s^2 + 36*J*K + 12*J*K*p + 12*K*k*p - 12*K*k*s + 12*K*k*p^2 - 12*K*k*s^2 + 15; %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tbl.Properties.VariableNames = {'scaleRatio','K','steps','flops','relFlopsBern','relFlopsLagr'};
    % Print to file:
    fileName = sprintf('%s_%s',fileRoot,strjoin(regexp(bases(1,j).getName,'\d+|^[A-Z]+|\((\w+)\)','match'),'_'));
    writetable(tbl,[fileName '.dat'],'Delimiter','\t');
%     % Export figure:
%     legend('Location','Best')
%     saveas(gcf,fileName,'png')
end