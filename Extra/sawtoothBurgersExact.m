function U = sawtoothBurgersExact(t,x,a)
% Evaluates the exact solution to the inviscid Burgers equation with
% periodic boundary conditions and a "triangle wave" initial condition:
%
%  u(0,x) = a + 1 - abs(x)
%
% This exact solution is valid for -inf < x < inf and t > 0 (i.e. can be
% sampled beyond its breaking time). Deduced via M.O.C., conservation
% principles, theorem 16.14 in Smoller, 1994 (i.e. entropy) and some
% "careful experimentation" to determine the damping of the shockwave - see
% Thomas, 1999 (page 106).
%
% Arguments
%  t: array of instants
%  x: array of positions
%  a: array of background advection speeds (default: 1/2)
%
% Miquel Herrera, 3/12/2020
%
% Defaults and checks:
if nargin < 3
    a = 0;
else
    validateattributes(a,{'numeric'},{'nonnegative','finite'})
end
% Preprocess samples:
[X,T,U] = meshgrid(x,t,a); % preallocate U and store A "in one fell swoop"
X = mod(X + 1 - (T < 1).*U.*T - (T >= 1).*(U + (U + .5).*(T - 1)),2) - 1;
% Sample at T < 1, flattening line:
s = T < 1 & T >= X;
U(s) = U(s) + (1 + X(s))./(1 + T(s));
% Idem, steepening line:
s = T < 1 & T < X;
U(s) = U(s) + (1 - X(s))./(1 - T(s));
% Sample remaining:
s = T >= 1;
U(s) =  .5 + U(s) + X(s)./(1 + T(s));
end