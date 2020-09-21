function f = objFun_peakPosition(basis,r)
% The smaller the dispersion over/undershoot peak, and the further to the 
% right it is located, the better (assigns lower cost).
%
% I/O
%  basis: discretization to evaluate
%  r: weighting factor; 0 -> peak to the right, 1 -> small over/undershoot
%  ---
%  f: cost <-> badness <-> penalty (smaller is better)
%
if nargin == 1
    r = .5; % peak strength and position are equally important
end
%
[z,k] = basis.getFourierFootprint;
w = 1i*z(1,k > 0);
k(k <= 0) = [];
[y,id] = max(real(w));
f = r*abs(y - k(id)) + (1-r)*(k(end) - k(id));
end