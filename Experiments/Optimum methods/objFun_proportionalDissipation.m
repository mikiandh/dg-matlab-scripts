function f = objFun_proportionalDissipation(basis,r,cutoff)
% Assigns a lower cost to the dispersion and dissipation relations that
% minimize the L2 norm of their difference, respectively, with: 
%  1) exact dispersion relation
%  2) target dissipation relation (proportional to the exact one for 
%     well-resolved wavenumbers, but to the dissipation error afterwards).
%
% I/O
%  basis: discretization to evaluate
%  r: ratio between dissipation and dispersion errors (after cutoff)
%  cutoff: use cutoff wavenumber? (default: true)
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
if nargin < 2
    r = 1; % r = 0 -> no dissipation; r >> 1 -> spectral cutoff filter
end
if nargin < 3
    cutoff = true;
end
[z,k] = basis.getFourierFootprint;
w = 1i*z(1,k > 0);
k(k <= 0) = [];
err = real(w) - k;
nc = find(err >= 0,1,'last'); % cutoff wavemode (right-most intercept)
if ~cutoff || isempty(nc)
    nc = 1; % add dissipation from the beginning
end
w0 = k - r*1i*(k >= k(nc)).*abs(err); % target Block wave
f = norm(w - w0);
end