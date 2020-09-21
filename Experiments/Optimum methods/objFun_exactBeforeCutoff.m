function f = objFun_exactBeforeCutoff(basis,cutoff)
% Assigns a lower cost to the dispersion and dissipation relations that
% minimize the L2 norm of their difference wi th their exact counterparts,
% but only in the range 0 < k < cutoff.
%
% I/O
%  basis: discretization to evaluate
%  cutoff: proportion of the spectrum that is targeted (fraction of 1)
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
if nargin < 2
    cutoff = .5; % half of the wavenumbers
end

[z,k] = basis.getFourierFootprint;
isTargeted = k > 0 & k < pi*basis.basisCount*cutoff;
f = norm(1i*z(1,isTargeted) - k(isTargeted));
end