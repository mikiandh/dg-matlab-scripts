function data = toro1(data,fileNameRoot)
% Toro's test problem 1 in p. 372, a transonic version of Sod's shock tube.

    function y = exactSolution(t,x)
        [r,u,p] = riemannEulerExact(t,x,1,0.75,1,0.125,0,0.1,0.3);
        y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
    end

% Preprocess:
norms = Norm({'ErrorL1','TV','BaselineTV'});
bcs = Farfield(exactSolution(0,0),exactSolution(0,1));
mesh = Mesh(data.basis,[0 1],bcs,data.K);
solver = SSP_RK4_10(Euler,[0 .2],...
    'norm',norms,...
    'exactSolution',@exactSolution,...
    'limiter',data.limiter);
solver.courantNumber = data.relCFL*solver.optimizeCFL(data.basis);

% Solve:
solver.initialize(mesh)
solver.launch(mesh)

% Postprocess:
data.wallClockTime = solver.wallClockTime;

data.densityErrorL1 = norms(1).vals(1);
data.momentumErrorL1 = norms(1).vals(2);
data.energyErrorL1 = norms(1).vals(3);

data.densityTV = norms(2).vals(1);
data.momentumTV = norms(2).vals(2);
data.energyTV = norms(2).vals(3);

data.densityExactTV = norms(3).vals(1);
data.momentumExactTV = norms(3).vals(2);
data.energyExactTV = norms(3).vals(3);

data.sensorRatio = solver.limiters(1).sensor.cumulativeActivationRatio;
data.limiterRatio = solver.limiters(1).cumulativeActivationRatio;

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],32)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end