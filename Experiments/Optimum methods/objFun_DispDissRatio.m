function f = objFun_DispDissRatio(basis)
% Objective function that penalizes the ratios between dispersion and
% dissipation the further they deviate from unity.
%
% I/O
%  basis: discretization to evaluate
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
% Mod. wavenumbers:
[z,k,e] = basis.getFourierFootprint;
z(:,k < 0) = [];
z = 1i*z;
e(:,k < 0) = [];
k(k < 0) = [];
% Ratios:
r = -imag(z)/basis.basisCount + 1e-3;
for i = 1:basis.basisCount
    r(i,:) = (abs(gradient(real(z(i,:)),k)-1) + 1e-3)./r(i,:);
end
% Energy-weighted L1 norm:
f = sum(trapz(k,abs(1 - r).*e,2))/basis.basisCount^2;
end