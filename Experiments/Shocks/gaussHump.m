function data = gaussHump(data,fileNameRoot,ptSkip)
% Analytic (i.e. as smooth as it gets) yet multichromatic initial condition
% for linear advection.

if nargin < 3
    ptSkip = 32;
end

physics = Advection(1);

    function y = exactSolution(t,x)
        y = physics.MOC(t,x,@(eta) exp(-9*pi/4*eta.^2),[-1 1]);
    end

% Preprocess:
norms = Norm({'ErrorL1','TV','BaselineTV'});
mesh = Mesh(data.basis,[-1 1],Periodic(2),data.K);
solver = SSP_RK4_10(physics,[0 8],...
    'norm',norms,...
    'exactSolution',@exactSolution,...
    'limiter',data.limiter);
if isnan(data.relCFL)
    solver.timeDelta = 1e-2;
    solver.isTimeDeltaFixed = true;
else
    solver.courantNumber = data.relCFL*solver.optimizeCFL(data.basis);
end

% Solve:
solver.initialize(mesh)
solver.launch(mesh)

% Postprocess:
data.wallClockTime = solver.wallClockTime;

data.densityErrorL1 = norms(1).vals(1);

data.densityTV = norms(2).vals(1);

data.densityExactTV = norms(3).vals(1);

data.sensorRatio = solver.limiters(1).sensor.cumulativeActivationRatio;
data.limiterRatio = solver.limiters(1).cumulativeActivationRatio;

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],ptSkip)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end