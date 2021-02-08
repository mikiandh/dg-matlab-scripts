function data = jiangShu(data,fileNameRoot)
% Initial condition for linear advection combining smooth, continuous and
% discontinuous profiles.
%
% Identical to that in Jiang and Shu (1996).
%
% Superposition of Gaussian and half-ellipse portions can be adjusted via
% weights.

physics = Advection(1);

a = 0.5;
z = -0.7;
delta = 0.005;
alpha = 10;
beta = log(2)/(36*delta^2);
wG = [1 4 1]/6;
wF = [1 4 1]/6;
    
    function y = exactSolution(t,x)
        y = physics.MOC(t,x,@initialCondition,[-1 1]);
    end

    function u = initialCondition(x)
        u = wG*[G(x,beta,z-delta); G(x,beta,z); G(x,beta,z+delta)].*inRange(x,-0.8,-0.6);
        u = u + 1*inRange(x,-0.4,-0.2);
        u = u + (1-abs(10*(x - 0.1))).*inRange(x,0,0.2);
        u = u + wF*[F(x,alpha,a-delta); F(x,alpha,a); F(x,alpha,a+delta)].*inRange(x,0.4,0.6);
    end

    function g = G(x,beta,z)
        g = exp(-beta*(x-z).^2);
    end

    function f = F(x,alpha,a)
        f = sqrt(max(1 - alpha^2*(x - a).^2,0));
    end

    function h = inRange(x,x1,x2)
        h = heaviside(x-x1) - heaviside(x-x2);
    end

% Preprocess:
norms = Norm({'ErrorL1','TV','BaselineTV'});
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

data.densityExactTV = norms(3).vals(1);

data.sensorRatio = solver.limiters(1).sensor.cumulativeActivationRatio;
data.limiterRatio = solver.limiters(1).cumulativeActivationRatio;

% Export:
solver.writeSolutionToFile([fileNameRoot '_solution'],32)
solver.writeLimiterToFile([fileNameRoot '_limiter'])
% savefig(sprintf('%s.fig',fileNameRoot))
end