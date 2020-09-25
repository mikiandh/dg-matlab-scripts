function f = objFun_Adams2015(basis)
% Objective function that tries to make dissipation proportional to
% dispersion; inspired by Adams, 2015.
%
% I/O
%  basis: discretization to evaluate
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
[z,k] = basis.getFourierFootprint;
z = 1i*z(1,k > 0);
k(k < 0) = [];
[~,kc] = resolvingEfficiency(z,k,basis.basisCount);
r = (abs(gradient(real(z),k) - 1)+1e-3)./(abs(imag(z))+1e-3); % current disp/diss ratios
r0 = max(r)*(k < kc); % target disp/diss ratios
f = trapz(abs(r - r0));
end