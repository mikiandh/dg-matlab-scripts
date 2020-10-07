function f = objFun_Asthana2015(basis)
% Objective function from Asthana & Jameson (2015).
%
% I/O
%  basis: discretization to evaluate
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
[z,k,e] = basis.getFourierFootprint('order','energy');
f = trapz(k,abs(1 - exp(1i*100*(k - 1i*z))).*e,2);
f = sum(f)/(2*basis.basisCount*pi); % I prefer this normalization
end