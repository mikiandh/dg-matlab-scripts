function [dPn, Pn] = diffLegendrePolynomialTable(n,x,Pn)
% Returns a matrix of evaluations of the first derivative of the Legendre 
% polynomials of degrees from 0 to n at x. Calculates it recursively. Uses
% a precomputed Pn if available, or computes it if not.
%
% Arguments
%  n: polynomial degree array
%  x: evaluation locations (1D array)
%
% Return
%  dPn: length(x) x n matrix
%
x = x(:); % make column vector
dPn = zeros(length(x),n); % preallocate
if nargin == 2 && n > 0
    Pn = legendrePolynomialTable(n,x); % compute necessary P_k(x)
end
% P'_0:
dPn(:,1) = 0;
% P'_k, k > 0:
for i = 1:n
    dPn(:,i+1) = i*Pn(:,i) + x.*dPn(:,i);
end
end