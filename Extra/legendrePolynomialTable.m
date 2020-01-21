function Pn = legendrePolynomialTable(n,x)
% Returns a matrix of evaluations of the Legendre polynomials of degrees
% from 0 to n at x. Calculates it recursively.
%
% Arguments
%  n: greatest polynomial degree to evaluate (n >= 0)
%  x: evaluation locations (1D array)
%
% Return
%  Pn: length(x) x n matrix
%
x = x(:); % make column vector
Pn = zeros(length(x),n+1); % preallocate
% P_0:
Pn(:,1) = 1;
if n == 0
    return;
end
% P_1:
Pn(:,2) = x;
% P_k, k > 1:
for i = 2:n
    Pn(:,i+1) = ((2*i-1)*x.*Pn(:,i) - (i-1)*Pn(:,i-1))/i;
end
end