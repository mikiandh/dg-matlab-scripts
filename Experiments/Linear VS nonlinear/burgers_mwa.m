clc
clear
close all

% This script computes a posteriori modified wavenumbers for the Burgers
% equation, using Stefan et al. (2014)'s approach.

%% Inputs
physics = Burgers;
basis = DGSEM(2);
K = 100; % number of patches
M = 10000; % number of samples, mesh-wide (a shitton of them)
L = 1; % signal length
isExport = false;

%% Safety checks
for aux = {K,M}
    validateattributes(aux{1},{'double'},{'even'})
end

%% Setup
solver = SSP_RK3(physics,[0 0]);
mesh = Mesh(basis,[-1 1]*L/2,Periodic(2),K);
t = linspace(mesh.edges([1 end]).coord,M+1); t(end) = []; % sample points
dx = mesh.elements(1).dx; % cell width
J = mesh.bases.basisCount; % DOFs/cell
N = K*J/2; % Nyquist wavemode

%% Wavemode loop
wavenums = zeros(3,N); % preallocation
D = parallel.pool.DataQueue;
D.afterEach(@plotMWA);
% Status:
fprintf('%s\n J = %d, K = %d\n',basis.getName,J,K)
fprintf(' Nyquist wavemode: %d\n',N)
fprintf(' Number of samples: %d\n\n',M)
tic
parfor n = 1:N
    % Define some useful quantities:
    wavelen = L/n; % wavelength
    wavenum = 2*pi/wavelen; % (unmodified) wavenumber
    fun = @(x) 1+.1*sin(2*pi*x/wavelen); % initial condition
    
    % Sample solution and residuals:
    solver.initialize(mesh,'initialCondition',fun), close %#ok<PFBNS>
    y0 = fun(t); % exact initial condition samples
    y = mesh.sample(t); % numerical initial condition samples
    z = mesh.sampleResidual(t); % numerical residual samples
    
    % Transform to Fourier domain:
    f = [0:M/2-1 -M/2:-1]; % wavenumbers
    fy0 = fft(y0); % exact Fourier initial condition
    fy = fft(y); % numerical Fourier initial condition
    fz0 = -1i*2*pi*f.*fy0; % exact Fourier residuals (spectral diff.)
    fz = fft(z); % numerical Fourier residuals
    
    % Reconstruct exact residuals:
    z0 = ifft(fz0);
    
    % Plot:
    send(D,{n,N,M,t,y0,y,z0,z,f,fy0,fy,fz0,fz})
    
    % Compute wavenumbers:
    wavenums(:,n) = [
        wavenum
        1i*fz0(n+1)./fy0(n+1)
        1i*fz(n+1)./fy(n+1)
        ];
end
toc

%% Dispersion and dissipation
aux = wavenums*dx/J; % scaled and dimensionless wavenumbers
figure(2)

% Dispersion:
subplot(2,1,1)
plot(aux(1,:),real(aux(2,:)),'DisplayName',sprintf('%s (%s)','Spectral method',solver.physics.getInfo))
hold on
plot(aux(1,:),real(aux(3,:)),'DisplayName',sprintf('%s (%s)',basis.getName,solver.physics.getInfo))
hold off
ylabel('$\tilde{\kappa}_{\Re} \Delta x / J$','Interpreter','latex')

% Dissipation:
subplot(2,1,2)
plot(aux(1,:),imag(aux(2,:)),'DisplayName',sprintf('%s (%s)','Spectral method',solver.physics.getInfo))
hold on
plot(aux(1,:),imag(aux(3,:)),'DisplayName',sprintf('%s (%s)',basis.getName,solver.physics.getInfo))
hold off
ylabel('$\tilde{\kappa}_{\Im} \Delta x / J$','Interpreter','latex')

% Finishing touches:
for i = 1:2
    subplot(2,1,i)
    xlabel('$\kappa \Delta x / J$','Interpreter','latex')
    legend('Location','Best')
end

%% Export
if ~isExport
    return
end
% Arrange data in a table:
tbl = array2table(aux.');
tbl.k = tbl{:,1};
tbl.k_real_spectral = real(tbl{:,2});
tbl.k_imag_spectral = imag(tbl{:,2});
tbl.k_real = real(tbl{:,3});
tbl.k_imag = imag(tbl{:,3});
tbl(:,1:3) = [];

% Save table to file:
fileName = sprintf('%s_%d_%s_%s','MWA',K,class(physics),strjoin(regexp(basis.getName,'\d+|^[A-Z]+|\((\w+)\)','match'),'_'));
writetable(tbl,[fileName '.dat'],'Delimiter','\t')

% Save also the figure:
savefig(fileName)