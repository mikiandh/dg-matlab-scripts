function data = woodwardColella(data,fileNameRoot)
% Woodward and Colella's (1984) interaction between two blast waves.

    function y = initialSolution(x)
        y = ones(3,length(x));
        y1 = Euler.primitiveToState([1 0 1000]');
        y2 = Euler.primitiveToState([1 0 0.01]');
        y3 = Euler.primitiveToState([1 0 100]');
        y(:,x < 0.1) = y1.*y(:,x < 0.1);
        y(:,x >= 0.1 & x < 0.9) = y2.*y(:,x >= 0.1 & x < 0.9);
        y(:,x > 0.9) = y3.*y(:,x > 0.9);
    end

% Preprocess:
data.K = max(data.K,200); % ensure a minimum of resolution
norms = Norm('TV');
mesh = Mesh(data.basis,[0 1],Reflective(0,0),data.K);
solver = SSP_RK4_10(Euler,[0 0.038],...
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