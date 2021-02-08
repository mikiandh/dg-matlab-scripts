function data = shuOsher(data,fileNameRoot)
% Shu-Osher's (1989) problem, mimicking a shockwave-turbulence interaction.

    function y = initialSolution(x)
        y = ones(3,length(x));
        y(:,x < -4) = [3.857143; 2.629369; 10.33333].*y(:,x < -4);
        y(1,x >= -4) = 1+0.2*sin(5*x(x >= -4));
        y(2,x >= -4) = 0.*y(2,x >= -4);
        y = Euler.primitiveToState(y);
    end

% Preprocess:
norms = Norm('TV');
mesh = Mesh(data.basis,[-5 5],Transmissive(2),data.K);
solver = SSP_RK4_10(Euler,[0 1.8],...
    'norm',norms,...
    'limiter',data.limiter);
solver.courantNumber = data.relCFL*solver.optimizeCFL(data.basis);

% Solve:
solver.initialize(mesh,'initialCondition',@initialSolution)
solver.launch(mesh)

% Postprocess:
data.wallClockTime = solver.wallClockTime;

data.densityTV = norms(1).vals(1);
data.momentumTV = norms(1).vals(2);
data.energyTV = norms(1).vals(3);

data.sensorRatio = solver.limiters(1).sensor.cumulativeActivationRatio;
data.limiterRatio = solver.limiters(1).cumulativeActivationRatio;

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],32)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end