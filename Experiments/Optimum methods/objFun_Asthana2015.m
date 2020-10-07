function f = objFun_Asthana2015(basis)
% Objective function a la Asthana & Jameson (2015).
%
% I/O
%  basis: discretization to evaluate
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
[z,k,e] = basis.getFourierFootprint;
z(:,k < 0) = [];
e(:,k < 0) = [];
k(k < 0) = [];
f = trapz(k,abs(1 - exp(1i*100*(k - 1i*z))).*e,2);
f = sum(f)/basis.basisCount^2;
end