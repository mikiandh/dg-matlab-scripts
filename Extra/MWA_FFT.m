function [kMod,k] = MWA_FFT(Nx,p,space,flag)
%% Modified wavenumber estimation
% All representable integer and positive wavenumbers in a FEM
% discretization.
%   kMin = 2*pi / L = 1 (L: domain size)
%   kMax = 2*pi / 2*l = #DOFs/2 (l = L/#DOFs: characteristic size)
%
if nargin < 4
    flag = false;
end
mesh = Mesh(linspace(0,2*pi,Nx+1),p,space);
k = 1:1:mesh.dofCount/2-1; % all positive wavenumbers
if flag
    fprintf(1,'%-8s\t %-8s\t %-8s\t %-8s\t %-8s\t %-8s\t %-8s\n',...
        'Method','Degree','#Patches','#Dimensions','#DOFs','k_min','k_max')
    fprintf(1,'%-8s\t %-8d\t %-8d\t %-8d\t %-8d\t %-8d\t %-8d\n',...
        class(space),p,Nx,mesh.bases.basisCount,mesh.dofCount,k(1),k(end))
end

% Sample locations: %%% TO DO: it seems NECESSARY to sample at quadrature points, but then FFT only works for p < 2 and NFFT doesn't work at all!
% x = mesh.getGaussLocations; % (possibly) non-uniformly, element-wise
% x = linspace(0,2*pi,mesh.dofCount); % uniformly, not element-wise
% x = (2*pi/(Nx*Ns*2*(p+1)).*(1:2:(Nx*Ns*2*(p+1)))); % p+1 points uniformly inside each span (will oversample IGA)
% x = (2*pi/(Nx*Ns*2*(2*p+2)).*(1:2:(Nx*Ns*2*(2*p+2)))); % 2p+2 points uniformly inside each span (oversampling)
% x = (2*pi/(Nx*Ns*2*(1)).*(1:2:(Nx*Ns*2*(1)))); % 1 point at the centroid of each span (undersampling)
x = (2*pi/(Nx*2*mesh.bases.basisCount).*(1:2:(Nx*2*mesh.bases.basisCount))); % N_b points uniformly inside each patch
if mod(length(x),2)
    error('#DOFs must be even. Try again.')
end
% Set up a relation between a single input and all admissible output wavenumbers:
i = @(j) j; % physical one; somewhat a la Hu et al, 1999
%     function n = i(j) % all admissible ones (both physical and parasite)
%         n0 = j; % presumed physical one
%         n1 = j+Nx:Nx:k(end);
%         n2 = j-Nx:-Nx:0;
%         n3 = -j+Nx:Nx:k(end);
%         n = horzcat(n0,n1,n2,n3);
%         n(n<0) = [];
%         n = n(1:p+1);
%     end
% Preallocate modified wavenumbers:
kMod = complex(1i*nan(length(i(1)),length(k)));
% Set up a non-uniform fast Fourier transform (Keiner et al. 2009):
plan = infft(x./(2*pi)-.5); % samples need to be rescaled to [-0.5, 0.5)
    function fHat = nfft(f)
        plan.f = f;
        infft_trafo(plan);
        fHat = plan.fcheck;
        %infft_trafo_direct(plan);
        %fHat = plan.fcheck_direct;
        fHat = fHat.'*length(fHat);
        %fHat(1) = fHat(1)'; % mode -N/2 becomes mode N/2
        fHat = ifftshift(fHat); % order is now that of MATLAB's FFT
    end
% Set up a non-equispaced discrete Fourier transform (by Vilnis Liepins):
%     function fHat = nfft(f)
%         N = length(x);
%         fn = [-ceil((N-1)/2):floor((N-1)/2)]/N; % NO IDEA; DOES NOT WORK
%         fHat = nedft(f,x,fn);
%         fHat = ifftshift(fHat);
%     end

% Wavenumber loop
%  Will measure the wavenumbers that 'appear' on the discretization
%  (output), when 'excited' with a monochromatic signal (input). Loops over
%  all integer wavenumbers resolvable by the discretization. Only 
%  positive-real-part modified wavenumbers are considered 
%  (single sided Fourier spectrum).

%%% IDEA: simply skip the wavenumbers where overlap occurs (thus avoiding the spikes)
% z = k(end)/mesh.bases.basisCount;
% ids = z:z:k(end);
% z = k;
% z(ids(1:end-1)) = [];
for j = k %z
%%%
    mesh.bases.project(mesh,[],@(x) sin(j*x))
    mesh.computeResiduals(Advection);
    % Sample states and residuals:
    [states,residuals] = mesh.sampleDOFs(x,'states','residuals');
    % Discrete Fourier transform:
    num = fft(residuals);
    den = fft(states);
    % Throw away negative wavenumbers:
    num(k(end)+2:end) = [];
    den(k(end)+2:end) = [];
    % Retrieve the modified wavenumber (1 among all the 'excited' ones):
    if i(j) == length(num)
        break % stop early if undersampling
    else
        kMod(:,j) = 1i*(num(i(j)+1)./den(i(j)+1));
    end
    % Retrieve modified wavenumbers: %%% NOT WORKING YET %%%
%     i = j + (-k(end):Nx:k(end)); % all those which SHOULD be excited (Van den Abeele, 2009)
%     i = i(i < k(end));
%     i = i(i > 0);
%     i = abs(num) > 1e-10 & abs(den) > 1e-10; % those which are detected as excited (FAILS)
%     kMod(:,j) = 1i*(num(i)./den(i));
    if flag
        y = linspace(0,2*pi,mesh.dofCount*30);
        [states_fine,residuals_fine] = mesh.sampleDOFs(y,'states','residuals');
        % Plot states:
        subplot(2,2,1)
        plot(y,sin(j*y),'k--')
        hold on
        plot(x,states,'bx')
        plot(y,states_fine,'b-')
        hold off
        xlabel('x')
        ylabel('Solution')
        xlim([0 2*pi])
        % Plot residuals:
        subplot(2,2,2)
        plot(x,residuals,'rx')
        xlabel('x')
        ylabel('Residual')
        hold on
        plot(y,residuals_fine,'r-')
        hold off
        xlim([0 2*pi])
        % Plot Fourier transform of states:
        den = den./k(end);
        subplot(2,2,3)
        bar([0 k],abs(den),'b')
        xlabel('^{|\Omega|}/_{\lambda}')
        ylabel('Magnitude of FFT of solution')
%         subplot(3,2,5)
%         plot([0 k],unwrap(angle(den)),'b:+')
%         xlabel('\xi')
%         ylabel('Angle of FFT of solution')
        % Plot Fourier transform of residuals:
        num = num./k(end);
        subplot(2,2,4)
        bar([0 k],abs(num),'r')
        xlabel('^{|\Omega|}/_{\lambda}')
        ylabel('FFT of residual')
%         subplot(3,2,6)
%         plot([0 k],unwrap(angle(num)),'r:+')
%         xlabel('\xi')
%         ylabel('Angle of FFT of residual')
        drawnow limitrate
        pause
    end
end
% Normalize from 0 to pi:
k = 2*pi*k/mesh.dofCount;
kMod = 2*pi*kMod/mesh.dofCount;
% 'Disguise' outliers (i.e. spikes):
% kMod = filloutliers(real(kMod),'linear','movmedian',Nx) +...
%     1i*filloutliers(imag(kMod),'linear','movmedian',Nx);
end