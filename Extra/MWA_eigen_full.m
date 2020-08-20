function [kMod,k] = MWA_eigen_full(Nx,p,disc,upwind)
%% Modified wavenumber estimation
% Compute modified wavenumbers via the eigendecomposition approach (van den
% Abeele, 2009, B1.1). Returns all modified wavenumbers (both physical and
% spurious) nondimensionalized by the patch size (i.e. not normalized).
%
if nargin < 4
    upwind = 1;
end
disc = disc.clone(p);
N = Nx*disc.basisCount; % #DOFs
N = round(N/2); % positive wavenumbers only
% Set up:
invMass = inv(disc.massMatrix);
% Modified wavenumber loop:
I = disc.basisCount;
k = 2*pi/Nx*(1:N); % all positive (nondimensional) wavenumbers
kMod = nan(I,N); % eigenvalues
for j = 1:N
    % Assemble discretization operator:
    if isa(disc,'FR')
        M = 2*disc.gradientMatrix'...
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
    % Retrieve modified wavenumbers:
    D = 1i*eigs(M,disc.basisCount); % 1i*eigenvalue <--> wavenumber
    % Apply a consistent re-ordering to eigenvalues:
    D = sortEigenvalues(D,j == 1);
    % Store them:
    % k(j) = k(j)/disc.basisCount; % normalize to [0,pi]
    % kMod(:,j) = D/disc.basisCount;
    kMod(:,j) = D;
end
end

%% Sorting routine
% Sort a column array of complex eigenvalues such that each one is closest to
% the previous instance of itself (stored as a persistent variable).
function [k,ids] = sortEigenvalues(k,reset)
% Set up the inital case:
persistent k0
if reset
    ids = length(k):-1:1; % strictly speaking, no order is enforced. In practice, this seems to always put the physical mode first.
    k = k(ids);
    k0 = k;
    return
end
% Determine a consistent order:
ids = 1:length(k);
for j = 1:length(k)
    dk = abs(k - k0(j));
    aux = min(dk);
    for i = 1:length(k)
        if dk(i) == aux
            ids(j) = i;
            break
        end
    end
end
% Apply it to the current set of eigenvalues:
k = k(ids);
k0 = k;
end