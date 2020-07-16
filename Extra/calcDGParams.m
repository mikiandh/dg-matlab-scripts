function [n,k,p] = calcDGParams(n,k,p)
% Given 2 out of 3 DG discretization parameters, this functions computes
% the remaining one. When a conflict occurs, the lowest value of n is
% preferred.
%
% Consistency check:
if sum(isnan([n k p])) > 1
    error('Discretization undefined (too many unknowns).')
end
% Case-by-case: 
if isnan(k)
    k = max(1,floor(n/(p+1)));
elseif isnan(p)
    p = max(0,floor(n/k-1));
end
% Correction:
n = k*(p+1);