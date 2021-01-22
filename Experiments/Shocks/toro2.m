function data = toro2(data,fileNameRoot)
% Toro's test problem 2 in p. 372, the 1-2-3 problem (symmetric expansion).

    function y = exactSolution(t,x)
        [r,u,p] = riemannEulerExact(t,x,1,-2,0.4,1,2,0.4,0.5);
        y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
    end

% Preprocess:
data.K = data.K + ~mod(data.K,2); % ensure 1 element is right at the middle
norms = Norm({'ErrorL1','TV'});
mesh = Mesh(data.basis,[0 1],Transmissive(2),data.K);
solver = SSP_RK3(Euler,[0 .15],...
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
data.momentumErrorL1 = norms(1).vals(2);
data.energyErrorL1 = norms(1).vals(3);

data.densityTV = norms(2).vals(1);
data.momentumTV = norms(2).vals(2);
data.energyTV = norms(2).vals(3);

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