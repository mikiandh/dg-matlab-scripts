function [kMod,k] = MWA_eigen(Nx,p,disc,upwind)
%% Modified wavenumber estimation
% Compute modified wavenumbers via the eigendecomposition approach (van den 
% Abeele, 2009, B1.1). Returns only the physical wavenumber, defined by
% Hu et al, 1999 as the one which is closest to the exact wavenumber over
% an extended range.
%
if nargin < 4
    upwind = 1;
end
disc = disc.clone(p);
N = Nx*disc.basisCount; % #DOFs
N = round(N/2); % positive wavenumbers only
% Set up
invMass = inv(disc.massMatrix);
% Modified wavenumber loop:
k = 2*pi/Nx*(1:N); % all positive (nondimensional) wavenumbers
kMod = nan(1,N); % eigenvalues
kOld = 0;
for j = 1:N
    % Assemble discretization operator:
    if isa(disc,'FR')
        M = 2*disc.gradientMatrix...
            + (1-upwind)*disc.correctionsR'*(-disc.right' + disc.left'*exp(1i*k(j)))...
            + (1+upwind)*disc.correctionsL'*(-disc.left' + disc.right'*exp(-1i*k(j)));
        M = -M;
    else
        M = 2*disc.gradientMatrix' +...
            (1-upwind)*disc.left*disc.left' -...
            (1+upwind)*disc.right*disc.right' +...
            exp(-1i*k(j))*(1+upwind)*disc.left*disc.right' -...
            exp(1i*k(j))*(1-upwind)*disc.right*disc.left';
        M = invMass*M;
    end
    % Retrieve all eigenvalues:
    D = 1i*eigs(M,disc.basisCount);
    % Determine the 'physical' one:
    [~,ids] = min(abs(D - kOld));
    kOld = D(ids);
    % Store it:
    kMod(j) = kOld;
end
% Normalize to [0,pi]:
k = k/disc.basisCount;
kMod = kMod/disc.basisCount;
end