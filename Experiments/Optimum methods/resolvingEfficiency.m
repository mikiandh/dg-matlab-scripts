function [e1,kc] = resolvingEfficiency(varargin)
% Computes the resolving efficiency of a basis or footprint (Lele, 1992).
switch nargin
    case 1 % basis
        [z,k] = varargin{:}.getFourierFootprint;
        J = varargin{:}.basisCount;
        z = 1i*z(1,k > 0);
        k(k < 0) = [];
    case 3
        [z,k,J] = varargin{:};
    otherwise
        error('Wrong inputs.')
end
n = find(imag(z) <= -0.01,1,'first'); % '1% rule' cutoff wavemode
kc = k(n);
e1 = kc/pi/J; % resolving efficiency
end