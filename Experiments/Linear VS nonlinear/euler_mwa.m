clc
clear
%close all

% This script computes a posteriori modified wavenumbers for the Euler
% equations, using Stefan et al. (2014)'s approach.

%% Inputs
Ma = .5; % Mach number
K = 100; % number of patches (i.e. resolution of disp. and diss. plots)
M = 10000; % number of samples, mesh-wide (a shitton, for good measure)
L = 1; % signal length
isExport = true; % export final plot and data?

for aux = {K,M}
    validateattributes(aux{1},{'double'},{'even'})
end

for basis = [
        DGSEM(2)
%         FR({'eta',0.0227795117943271},2)
%         DGIGA(1,2,1)
%         DGIGA_nodal(1,2,1)
%         DGSEM(5)
%         FR({'eta',0.0580283806157087},5)
%         DGIGA(2,3,1)
%         DGIGA_nodal(2,3,1)
%         DGSEM(19)
%         FR({'eta',0.152203674556189},19)
%         DGIGA(2,14,9)
%         DGIGA_nodal(2,14,9)
        ]'
    for physics = [
            Euler('Exact')
%             Euler('LLF')
%             Euler('RoeNoFix')
%             Euler('RoeHartenHyman')
%             Euler('HLL')
%             Euler('HLLE')
%             Euler('HLLC')
%             Euler('KEP')
            ]'
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
        fprintf(1,'%s\n J = %d, K = %d\n',basis.getName,J,K)
        fprintf(1,' Nyquist wavemode: %d\n',N)
        fprintf(1,' Number of samples: %d\n',M)
        tic
        parfor n = 1:N
            % Define some useful quantities:
            wavelen = L/n; % wavelength
            wavenum = 2*pi/wavelen; % (unmodified) wavenumber
            %%%%%%%% Initial condition %%%%%%%%%%%%%%%%
            primitiveFun = @(x) [
                1+.1*sin(2*pi*x/wavelen)
                ones(size(x))
                ones(size(x))/(1.4*Ma^2)
                ];
            stateFun = @(x) [
                1+.1*sin(2*pi*x/wavelen)
                1+.1*sin(2*pi*x/wavelen)
                .5*(1+.1*sin(2*pi*x/wavelen))+1/(.56*Ma^2)
                ];
            fluxFun = @(x) [
                1+.1*sin(2*pi*x/wavelen)
                1+.1*sin(2*pi*x/wavelen) + 1/(1.4*Ma^2)
                .5*(1+.1*sin(2*pi*x/wavelen))+1/(.56*Ma^2) + 1/(1.4*Ma^2)
                ];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sample solution and residuals:
            solver.initialize(mesh,'initialCondition',stateFun) %#ok<PFBNS>
            close
            y0 = stateFun(t); % exact initial condition samples
            y = mesh.sample(t); % numerical initial condition samples
            z0 = fluxFun(t); % exact fluxes
            z = mesh.sampleResidual(t); % numerical residual samples
            
            % Transform to Fourier domain:
            f = [0:M/2-1 -M/2:-1]; % wavenumbers
            fy0 = fft(y0,[],2); % exact Fourier initial condition
            fy = fft(y,[],2); % numerical Fourier initial condition
            fz0 = fft(z0,[],2); % exact Fourier fluxes
            fz0 = -1i*2*pi/L*f.*fz0; % exact Fourier residuals (spectral diff.)
            fz = fft(z,[],2); % numerical Fourier residuals
            
            % Reconstruct exact residuals:
            z0 = ifft(fz0,[],2,'symmetric');
            
            % Plot:
            send(D,{n,N,M,t,y0(1,:),y(1,:),z0(1,:),z(1,:),f,fy0(1,:),fy(1,:),fz0(1,:),fz(1,:)})
            
            % Compute wavenumbers:
            wavenums(:,n) = [
                wavenum
                1i*fz0(1,n+1)./fy0(1,n+1)
                1i*fz(1,n+1)./fy(1,n+1)
                ];
        end
        fprintf(1,'Finished in %g seconds.\n\n',toc)
        
        %% Dispersion and dissipation
        aux = wavenums*dx/J; % scaled and dimensionless wavenumbers
        figure(2)
        
        % Dispersion:
        subplot(2,1,1)
        hold on
        plot(aux(1,:),real(aux(2,:)),'DisplayName',sprintf('%s (%s)','Spectral method',solver.physics.getInfo))
        plot(aux(1,:),real(aux(3,:)),'DisplayName',sprintf('%s (%s)',basis.getName,solver.physics.getInfo))
        hold off
        ylabel('$\tilde{\kappa}_{\Re} \Delta x / J$','Interpreter','latex')
        
        % Dissipation:
        subplot(2,1,2)
        hold on
        plot(aux(1,:),imag(aux(2,:)),'DisplayName',sprintf('%s (%s)','Spectral method',solver.physics.getInfo))
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
            continue
        end
        % Arrange data in a table:
        % tbl = array2table(wavenums.');
        tbl = array2table(aux.');
        tbl.k = tbl{:,1};
        tbl.k_real_spectral = real(tbl{:,2});
        tbl.k_imag_spectral = imag(tbl{:,2});
        tbl.k_real = real(tbl{:,3});
        tbl.k_imag = imag(tbl{:,3});
        tbl(:,1:3) = [];
        
        % Save table to file:
        fileName = sprintf('%s_K=%d_J=%d_%s_%s',...
            'MWA',K,basis.basisCount,physics.getInfo,class(basis));
        fileName = regexprep(fileName,{'\s','(',')'},{'_','',''});
        writetable(tbl,[fileName '.dat'],'Delimiter','\t')
        
        % Save also the figure:
        savefig(fileName)
    end
end