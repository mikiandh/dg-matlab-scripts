function [kMod,k] = MWA_DSDFT(Nx,p,space,q,flag)
%% Modified wavenumber estimation
% All representable integer and positive wavenumbers in a FEM
% discretization.
%   kMin = 2*pi / L = 1 (L: domain size = 2*pi)
%   kMax = 2*pi / 2*l = #DOFs/2 (l = L/#DOFs: characteristic size)
%
if nargin < 4
    q = p+1; % quadrature degree (for the Fourier projection)
end
if nargin < 5
    flag = false;
end
mesh = Mesh(linspace(0,2*pi,Nx+1),p,space);
k = 1:1:mesh.dofCount/2;
if flag
    fprintf(1,'%-8s\t %-8s\t %-8s\t %-8s\t %-8s\t %-8s\t %-8s\n',...
        'Method','Degree','#Patches','#Dimensions','#DOFs','k_min','k_max')
    fprintf(1,'%-8s\t %-8d\t %-8d\t %-8d\t %-8d\t %-8d\t %-8d\n',...
        class(space),p,Nx,mesh.bases.basisCount,mesh.dofCount,k(1),k(end))
end
% Wavenumber loop
%  For each input wavenumber (j), the associated monochromatic signal is
%  projected onto the discretization, and its residuals are evaluated. The
%  approximate solution (piecewise polynomial in some basis) is then
%  projected onto a Fourier basis with a single wavenumber (i = +j).
%
i = @(j) j; % physical one; somewhat a la Hu et al, 1999
%     function n = i(j) % all admissible ones (both physical and parasite)
%         n0 = j; % presumed physical one
%         n1 = j+Nx:Nx:k(end);
%         n2 = j-Nx:-Nx:0;
%         n3 = -j+Nx:Nx:k(end);
%         n = horzcat(n0,n1,n2,n3);
%         n(n<0) = [];
%         n = n(1:p+1)';
%     end
kMod = complex(1i*nan(length(i(1)),length(k))); % preallocation of modified wavenumbers
[gaussCoords,gaussWeights] = Legendre.quadratureGaussLegendre(q);
gaussCoords = gaussCoords';
for j = k
    mesh.bases.project(mesh,[],@(x) sin(j*x))
    mesh.computeResiduals(Advection);
    % Projection to a sparse Fourier basis:
    x = zeros(1,Nx*(q+1));
    states = x;
    residuals = x;
    num = 0;
    den = 0;
    id = 1:q+1;
    for element = mesh.elements
        if element.basis.isHybrid, error('DGIGA is not supported.'), end
        phiL = element.basis.sampleAt(gaussCoords);
        statesk = element.states*phiL;
        residualsk = element.residuals*phiL;
        xk = element.mapFromReference(gaussCoords);
        phiF = exp(-1i*i(j).*xk);
        num = num + (residualsk.*phiF)*gaussWeights;
        den = den + (statesk.*phiF)*gaussWeights;
        %%%
        if flag
            x(id) = xk;
            states(id) = statesk;
            residuals(id) = residualsk;
            id = id + q+1;
        end
        %%%
    end
    kMod(:,j) = 1i*num./den;
    % Generate plots, if requested:
    if flag
        aux = num/(2*pi);
        num = zeros(1,k(end)+1);
        num(i(j)+1) = aux;
        aux = den/(2*pi);
        den = zeros(1,k(end)+1);
        den(i(j)+1) = aux;
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
        subplot(2,2,3)
        bar([0 k],abs(den),'b')
        xlabel('^{|\Omega|}/_{\lambda}')
        ylabel('Magnitude of FFT of solution')
%         subplot(3,2,5)
%         plot([0 k],unwrap(angle(den)),'b:+')
%         xlabel('\xi')
%         ylabel('Angle of FFT of solution')
        % Plot Fourier transform of residuals:
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
end