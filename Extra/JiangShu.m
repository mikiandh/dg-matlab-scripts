function u = JiangShu(x,a,z,delta,alpha,beta)
% Function typically used as initial condition for benchmarking the
% performance of limiters in 1D (Jiang and Shu, 1996). Nomenclature 
% based on Krivodonova, 2007.
%
% Arguments
%  x: evaluation location (in reference space such that -1 < x < 1)
%  parameters: a, z, delta, alpha, beta (optional)

switch nargin
    case 1
        a = 0.5;
        z = -0.7;
        delta = 0.005;
        alpha = 10;
        beta = 770.163533955495;
    case 6
        % ok
    otherwise
        error('Incompatible arguments.')
end
u = 0.*x;
u = u + 1/6*(G(x,beta,z-delta) + G(x,beta,z+delta) + 4*G(x,beta,z)).*...
    inRange(x,-0.8,-0.6);
u = u + 1*inRange(x,-0.4,-0.2);
u = u + (1-abs(10*(x - 0.1))).*inRange(x,0,0.2);
u = u + 1/6*(F(x,alpha,a-delta) + F(x,alpha,a+delta) + 4*F(x,alpha,a)).*...
    inRange(x,0.4,0.6);
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