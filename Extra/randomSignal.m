function u = randomSignal(x,n,m,L,seed)
% Random singal made up of n smooth tones of random amplitude and frequency
% and m discontinuous jumps randomly distributed and with random amplitude.
%
% Arguments
%  x: evaluation location (in reference space such that -1 < x < 1)
%  n: number of sinusoidal components
%  m: number of heaviside components
% seed: seed for the random number generator

switch nargin
    case 1
        n = 3;
        m = 3;
        L = [x(1) x(end)];
    case 2 % smooth signal
        m = 0;
        L = [x(1) x(end)];
    case 3
        L = [x(1) x(end)];
    case 4
        % ok
    case 5
        rng(seed);
    otherwise
        error('Wrong number of arguments.')
end
% Smooth components:
y = rand(n,3);
u = y(:,1)'*sin(2*pi/diff(L)*y(:,2)*x + y(:,3));
% Sharp components:
y = rand(m,2);
y(:,2) = diff(L)*(y(:,2) - min(y(:,2)))/(max(y(:,2)) - min(y(:,2))) + L(1);
y(y(:,2) == 1 | y(:,2) == 0) = 0;
u = y(:,1)'*heaviside(x - y(:,2)) + u;
u = 2*(u - min(u))/(max(u) - min(u)) - 1;
end