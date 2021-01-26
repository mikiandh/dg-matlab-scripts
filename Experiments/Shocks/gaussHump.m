function data = gaussHump(data,fileNameRoot)
% Analytic (i.e. as smooth as it gets) yet multichromatic initial condition
% for linear advection.

physics = Advection(1);

    function y = exactSolution(t,x)
        y = physics.MOC(t,x,@(eta) exp(-9*pi/4*eta.^2),[-1 1]);
    end

% Preprocess:
data.limiter(:,2:end) = [];
norms = Norm({'ErrorL1','TV'});
mesh = Mesh(data.basis,[-1 1],Periodic(2),data.K);
solver = SSP_RK4_10(physics,[0 8],...
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

data.densityTV = norms(2).vals(1);

data.troubledDofs = 0;
data.limitedDofs = 0;
for element = mesh.elements
    data.troubledDofs =...
        data.troubledDofs + element.isTroubled(:,:,1)*numel(element.isLimited(:,:,1));
    data.limitedDofs =...
        data.limitedDofs + sum(sum(element.isLimited(:,:,1)));
end

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],32)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end