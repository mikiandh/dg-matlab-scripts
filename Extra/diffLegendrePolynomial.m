function dP = diffLegendrePolynomial(n,x,sym_flag)
% Returns the 1st derivative of Legendre polynomial of degree n, evaluated 
% at points x. Uses the recusion formula inside the domain, and checks for
% values at endpoints (1, -1); alternatively, uses symbolic
% differentiation.
%
% Inspired by:
% Bücker and Willkom, "First-order Derivatives of Associated Legendre
% Functions", Oxford, 2016.
%
% Inputs
%  n: Legendre polynomial degree
%  x: evaluation points
%  sym_flag (false): use symbolic math instead
%
% Check for optional argument:
if nargin == 2 || sym_flag == 0
    % Preallocate:
    dP = zeros(size(x));
    % Exclude edge points (-1 and 1):
    ids = find(x~=1 & x~=-1); % ids of locations not on the boundaries
    xIn = x(ids);
    % Interior points:
    dP(ids) = (n+1)./(1-xIn.^2).*(xIn.*legendrePolynomial(n,xIn)-legendrePolynomial(n+1,xIn));
    % Right boundary:
    dP_right = 0.5*n*(n+1);
    dP(x == 1) = dP_right;
    % Left boundary:
    dP_left = (-1)^(n+1)*n*(n+1)*0.5;
    dP(x == -1) = dP_left;
else
    syms y
    P = legendreP(n,y);
    dP = diff(P);
    dP = double(subs(dP,y,x));
end
end