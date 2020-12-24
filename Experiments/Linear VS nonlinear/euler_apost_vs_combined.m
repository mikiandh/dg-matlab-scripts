clc
clear
%close all

% This script compares the results of a posteriori MWA with combined-mode
% analysis.

%% Inputs
Ma = .5; % Mach number
K = 100; % number of patches (i.e. resolution of disp. and diss. plots)
M = 10000; % number of samples, mesh-wide (a shitton, for good measure)
L = 1; % signal length
T = [.1 1 10]'; % sample instants
isExport = true; % export final plot and data?

for aux = {K,M}
    validateattributes(aux{1},{'double'},{'even'})
end

for basis = [
        DGSEM(0)
        DGSEM(1)
        DGSEM(2)
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
            send(D,{n,N,M,K,t,y0(1,:),y(1,:),z0(1,:),z(1,:),f,fy0(1,:),fy(1,:),fz0(1,:),fz(1,:)})
            
            % Compute wavenumbers:
            wavenums(:,n) = [
                wavenum
                1i*fz0(1,n+1)./fy0(1,n+1)
                1i*fz(1,n+1)./fy(1,n+1)
                ];
        end
        fprintf(1,'Finished in %g seconds.\n\n',toc)
        
        %% Dispersion and dissipation errors
        aux = dx*wavenums;
        data.k_aprio = aux(1,:)/basis.basisCount;
        data.ang_aprio = (aux(1,:) - real(aux(3,:))).*T;
        data.amp_aprio = exp(imag(aux(3,:)).*T);
        [data.k_combi,data.ang_combi,data.amp_combi] = basis.getAngAmp(T);
        data.k_combi = data.k_combi/basis.basisCount;
        data.ang_combi(:,data.k_combi < 0) = [];
        data.amp_combi(:,data.k_combi < 0) = [];
        data.k_combi(data.k_combi < 0) = [];
        figure(2)
        % Dispersion error:
        subplot(2,1,1)
        plot(data.k_aprio,data.ang_aprio,'.')
        hold on
        plot(data.k_combi,data.ang_combi)
        hold off
        ylabel('$\Delta\Psi$','Interpreter','latex')
        title(basis.getName)
        
        % Dissipation factor:
        subplot(2,1,2)
        plot(data.k_aprio,data.amp_aprio,'.')
        hold on
        plot(data.k_combi,data.amp_combi)
        hold off
        ylabel('$1 - G$','Interpreter','latex')
        
        % Finishing touches:
        legNames = compose('%s, T = %g',["a posteriori","combined-mode"]',T')';
        for i = 1:2
            subplot(2,1,i)
            xlabel('$\kappa \Delta x / J$','Interpreter','latex')
            legend(legNames{:},'Location','Best')
        end
        
        %% Export
        if ~isExport
            continue
        end
        fileName = sprintf('%s_K=%d_J=%d_%s_%s',...
            'mwa_vs_combi',K,basis.basisCount,physics.getInfo,class(basis));
        fileName = regexprep(fileName,{'\s','(',')'},{'_','',''});
        savefig(fileName)
    end
end