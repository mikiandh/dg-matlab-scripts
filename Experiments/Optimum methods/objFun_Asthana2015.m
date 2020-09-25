function f = objFun_Asthana2015(basis)
% Objective function a la Asthana & Jameson (2015).
%
% I/O
%  basis: discretization to evaluate
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
[z,k] = basis.getFourierFootprint;
z = 1i*z(1,k > 0);
k(k < 0) = [];
f = trapz(abs(1 - exp(1i*100*(k - z))))/basis.basisCount;
end