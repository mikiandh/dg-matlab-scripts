function data = gaussHump(data,fileNameRoot)
% Analytic (i.e. as smooth as it gets) yet multichromatic initial condition
% for linear advection.

physics = Advection(1);

    function y = exactSolution(t,x)
        y = physics.MOC(t,x,@(eta) exp(-9*pi/4*eta.^2),[-1 1]);
    end

% Preprocess:
norms = Norm({'ErrorL1','TV'});
mesh = Mesh(data.basis,[-1 1],Periodic(2),data.K);
solver = SSP_RK3(physics,[0 8],...
    'norm',norms,...
    'exactSolution',@exactSolution,...
    'limiter',data.limiter{1});
solver.courantNumber = data.relCFL*solver.optimizeCFL(data.basis);

% Solve:
solver.initialize(mesh)
solver.launch(mesh)

% Postprocess:
data.wallClockTime = solver.wallClockTime;

data.densityErrorL1 = norms(1).vals(1);

data.densityTV = norms(2).vals(1);

data.troubledDofs = 0;
data.limitedDofs = 0;
for element = mesh.elements
    data.troubledDofs =...
        data.troubledDofs + element.isTroubled*numel(element.isLimited(:));
    data.limitedDofs =...
        data.limitedDofs + sum(element.isLimited(:));
end

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],8)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
end