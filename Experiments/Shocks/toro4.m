function data = toro4(data,fileNameRoot)
% Toro's test problem 4 in p. 372 (also problem 5 in p. 129), a frontal
% collision between two blast waves.
% It reproduces the shock-shock interaction in Woodward and Colella's 
% (1984) test problem.

    function y = exactSolution(t,x)
        [r,u,p] = riemannEulerExact(t,x,5.99924,19.5975,460.894,5.99242,-6.19633,46.095,0.4);
        y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
    end

% Preprocess:
norms = Norm({'ErrorL1','TV','BaselineTV'});
mesh = Mesh(data.basis,[0 1],Transmissive(2),data.K);
solver = SSP_RK4_10(Euler,[0 0.035],...
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