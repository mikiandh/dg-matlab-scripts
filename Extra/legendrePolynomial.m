function Pn = legendrePolynomial(n,x)
% Overrides MATLAB's built-in. Hopefully is faster by avoiding symbolic
% math. Vectorized.
%
% Arguments
%  n: scalar polynomial degree
%  x: evaluation locations (1D array)
%
% Return
%  Pn: column vector of samples at x of Legendre polynomial of degree n
%
if length(n) ~= 1
    Pn = legendreP(n,x); % TODO: implement a separate recursive formula for such cases
    return
end
k = 0:n;
x = (0.5*(x(:)-1)).^(k); % x matrix
binomials = @(k) nchoosek(n,k)*nchoosek(n+k,n);
B = arrayfun(binomials,k); % vector of product of binomials
A = x.*B; % auxiliary matrix (pre-summation)
Pn = sum(A,2);
end