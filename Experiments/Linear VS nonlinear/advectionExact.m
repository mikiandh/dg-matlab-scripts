function y = advectionExact(t,x,fun,a,L)
% Evolves an initial condition according to the advection equation, with
% periodic boundary conditions.
%
% Arguments
%  t: instant
%  x: position
%  fun: initial condition (at t = 0)
%  a: advection speed
%  L: domain boundary positions (2-element array)
%
x = x-t*a; % sample coordinates in "extended domain"
x = L(1) + mod(x-L(1),L(end)-L(1)); % corresponding coordinates in the original domain
y = fun(x); % sample values in "extended domain"
end