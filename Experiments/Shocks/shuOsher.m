function data = shuOsher(data,fileNameRoot)
% Shu-Osher's (1989) problem, mimicking a shockwave-turbulence interaction.

    function y = initialSolution(x)
        y = ones(3,length(x));
        for i = 1:length(x)
            if x(i) < -4
                y(:,i) = Euler.primitiveToState([3.857143 2.629369 10.33333]');
            else
                y(:,i) = Euler.primitiveToState([1+0.2*sin(5*x(i)) 0 1]');
            end
        end
    end

% Preprocess:
data.K = max(data.K,300); % ensure a minimum of resolution
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