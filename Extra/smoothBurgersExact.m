function u = smoothBurgersExact(t,x,fun0)
% Solves the inviscid Burgers' equation assumming smooth initial conditions
% and solution until the end time. Miquel Herrera, 18-02-2020.
%
% Arguments:
%  t: sample time instants
%  x: sample locations
%  u0: initial solution function
%
% Return:
%  u: exact solution, at time t (column: position)
%
% Ensure proper arrangement:
t = reshape(t,[],1);
x = reshape(x,1,[]);
% Set up:
u0 = repmat(fun0(x),length(t),1);
err = inf;
tol = 1e-6;
iter = 0;
maxIter = 1e2;
% Naive iteration:
while err > tol &&  iter < maxIter % check convergence
    iter = iter + 1;
    u = fun0(x - t.*u0); % update
    err = max(abs(u - u0)); % max norm error
    u0 = u;
end