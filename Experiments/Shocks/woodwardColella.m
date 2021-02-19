function data = woodwardColella(data,fileNameRoot,ptSkip)
% Woodward and Colella's (1984) interaction between two blast waves.

if nargin < 3
    ptSkip = 32;
end

    function y = initialSolution(x)
        y = ones(3,length(x));
        y1 = Euler.primitiveToState([1 0 1000]');
        y2 = Euler.primitiveToState([1 0 0.01]');
        y3 = Euler.primitiveToState([1 0 100]');
        y(:,x < 0.1) = y1.*y(:,x < 0.1);
        y(:,x >= 0.1 & x < 0.9) = y2.*y(:,x >= 0.1 & x < 0.9);
        y(:,x >= 0.9) = y3.*y(:,x >= 0.9);
    end

% Preprocess:
norms = Norm('TV');
mesh = Mesh(data.basis,[0 1],Reflective(0,0),data.K);
solver = SSP_RK4_10(Euler,[0 0.038],...
    'norm',norms,...
    'limiter',data.limiter);
mtd = 'project';
if isnan(data.relCFL)
    solver.timeDelta = 1e-4;
    solver.isTimeDeltaFixed = true;
%     mtd = 'interpolate';
else
    solver.courantNumber = data.relCFL*solver.optimizeCFL(data.basis);
end

% Solve:
solver.initialize(mesh,'initialCondition',@initialSolution,'method',mtd)
solver.launch(mesh)

% Postprocess:
data.wallClockTime = solver.wallClockTime;

data.densityTV = norms(1).vals(1);
data.momentumTV = norms(1).vals(2);
data.energyTV = norms(1).vals(3);

data.sensorRatio = solver.limiters(1).sensor.cumulativeActivationRatio;
data.limiterRatio = solver.limiters(1).cumulativeActivationRatio;

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],ptSkip)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end