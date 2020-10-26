classdef Basis < matlab.mixin.SetGet & matlab.mixin.Heterogeneous
    properties (Constant,Abstract)
        isNodal
        isModal
        isHybrid
    end
    properties
        degree
        order
        basisCount
        left
        right
        massMatrix % inner products between every basis function (row) and test function (column)
        gradientMatrix % inner products between every basis function (row) and test function first derivative (column)
        dofCoords % positions associated with each degree of freedom
        breakCoords % positions where the approximate solution experiences a reduction in smoothness
        nodeCoords
        controlCoords % abscissae of all control points in a patch (in reference patch coordinates)
        gaussCoords
        gaussWeights
        pairs % 2D array of test function/basis function pairs within mutual nonzero support; 1st row: test function, r; 2nd row: basis function, j; 3rd row: linear index; 4th row: transposed linear index
        edges % idem, but excludes pairs of repeated control points, i.e. j == r ("virtual" edges a la Kuzmin et al. 2012)
    end
    methods (Abstract)
        clone(prototype,degree) % required for the prototype instantiation pattern
        sampleAt(this,x) % sample basis components at given locations
        interpolate(this,element,fun0) % sample-based projection (non-conservative)
    end
    methods (Abstract, Access = {?Basis,?Element})
        computeResiduals(this,element,physics) % evaluate discrete spatial operator, i.e residual (a la Wang et al. 2013)
        getLegendre(this,element,j,i) % get chosen Legendre coefficients from an element (assumed to employ this basis)
        setLegendre(this,element,modes,i) % set all Legendre coefficients of an element (assumed to employ this basis)
        getLagrange(this,element,j,i) % get chosen Lagrange coefficients from an element (assumed to employ this basis)
        setLagrange(this,element,nodes,i) % set all Lagrange coefficients of an element (assumed to employ this basis)
    end
    methods
        %% L2 projection (conservative)
        function project(this,element,fun0)
            % Mass-preserving projection onto this basis. Inner products 
            % between each test function and the exact function are 
            % approximated using adaptive quadrature. All basis subclasses
            % share this general implementation (break span based).
            %
            % IMPORTANT: 'this' it is assumed to be the basis of both trial
            % and test spaces. For non-standard Galerkin (e.g. FR),
            % override or extend as necessary.
            %
            % Set up:
            I = length(fun0(0));
            N = this.basisCount;
            subs = [this.breakCoords(1:end-1); this.breakCoords(2:end)];
            % Vectorize right-hand-side integrand:
            fv = @(xi) repmat(fun0(element.mapFromReference(xi)),N,1).*repelem(this.sampleAt(xi),I,1);
            % Approximate integration:
            rhs = Algorithms.quadvgk(fv,subs,I*N);
            % Re-arrange and solve:
            element.states = reshape(rhs,I,N)/this.massMatrix;
        end
        %% DOF sparsity graph
        function computeSparsityGraph(this)
            % Fills the properties related to the sparsity graph of the
            % discretization (i.e. nonzero entries in the mass matrix).
            [r,j] = find(this.massMatrix);
            ind = sub2ind(this.basisCount*[1 1],r,j);
            dni = sub2ind(this.basisCount*[1 1],j,r);
            this.pairs = [r'; j'; ind'; dni'];
            this.edges = this.pairs(:,r~=j);
            %%% this.edges = [this.edges; 1:size(this.edges,2)];
        end
        %% Quadrature coordinates
        function x = getGaussCoords(this)
            % Retrieves quadrature point locations in reference patch 
            % coordinates. Default implementation (override as necessary).
            x = this.gaussCoords;
        end
        %% Oneliner info
        function name = getName(this)
            % Returns the name of this basis and its degree.
            name = sprintf('%s, p = %d',strrep(class(this),'_','-'),this.degree);
        end
        %% Verbose info
        function getInfo(this,varargin)
            % Prints information about the discretization.
            % Set up name-value input argument list:
            opts = inputParser;
            addParameter(opts,'generic',true,@islogical) % display generic info about the discretization
            addParameter(opts,'resolution',false,@islogical) % display its resolution limits (spectral space)
            addOptional(opts,'mode',nan,@isnumeric) % display the 'position' of the slected mode in spectral space
            opts.KeepUnmatched = true;
            parse(opts,varargin{:});
            % Print some generic information:
            if opts.Results.generic
                if this.isNodal
                    disp('Basis is nodal.')
                else
                    disp('Basis is not nodal.')
                end
                if this.isModal
                    disp('Basis is modal.')
                else
                    disp('Basis is not modal.')
                end
                if this.isHybrid
                    fprintf('Basis has both DG (inter-patch) and CG (intra-patch) polynomial coupling. \n')
                    if isa(this,'Bspline')
                        contClass = this.smoothness;
                        if this.nonzeroSpanCount == 1
                            contClass = inf;
                        elseif isinf(contClass) || isnan(contClass)
                            contClass = this.degree - 1;
                        end
                        disp('Basis and test function spaces are identical (Bubnov-Galerkin). Basis functions are B-splines. Solution is approximated as a spline curve.')
                        fprintf('> The (open) knot vector is: [%g',this.knots(1)), fprintf(' %g', this.knots(2:end)), fprintf(']\n')
                        fprintf('> There are %d knots (%d of which are breakpoints) and %d control points patch-wide.\n',this.knotCount,length(this.breakCoords),length(this.controlCoords))
                        fprintf('> Approximate solutions within the patch are made of %d piecewise polynomials of degree %d, and belong to the C^%d smoothness class.\n',this.nonzeroSpanCount,this.degree,contClass)
                        %fprintf('Basis has %d dimensions. Patch-wide approximations have %d degrees of freedom and employ %d quadrature points. \n',this.basisCount,this.basisCount,length(this.gaussCoords))
                    end
                else
                    disp('Basis has DG coupling with adjacent patches.')
                    if isa(this,'Legendre')
                        fprintf('Basis functions are Legendre polynomials of degree %d. Solution is approximated as a polynomial of degree %d. ',this.degree,this.degree)
                    elseif isa(this,'Lagrange')
                        fprintf('Basis functions are Lagrange polynomials of degree %d, centered at the nodes. Solution is approximated as a Lagrange interpolant of degree %d. ',this.degree,this.degree)
                    end
                    if isa(this,'FR')
                        fprintf('\n')
                        fprintf('> Test functions are Dirac delta functions, centered at the nodes (collocation).\n')
                        if isnumeric(this.param)
                            fprintf('> Correction functions are polynomials of degree %d, obtained by c = %g.\n',this.order,this.param)
                        else
                            fprintf('> Correction functions are polynomials of degree %d, corresponding to the %s variant of FR/CPR.\n',this.order,this.param)
                        end
                    else
                        fprintf('Basis and test function spaces are identical (Bubnov-Galerkin).\n')
                    end
                end
                fprintf('Basis has %d dimensions. Patch-wide approximations have %d degrees of freedom and employ %d (Gauss-Legendre) quadrature points. \n',this.basisCount,this.basisCount,length(this.gaussCoords))
            end
            % Print information about resolution limits:
            if opts.Results.resolution
                SCALES = {
                    'Wave mode (-)', 0, '1', [num2str(this.basisCount/2),'Nx'];
                    'Wave length (m)', Inf, 'Nx*dx', [num2str(2/this.basisCount),'dx'];
                    'Wave number (1/m)', 0, [num2str(2*pi),'/(Nx*dx)'], [num2str(this.basisCount*pi),'/dx'];
                    'Phase angle (º)', 0, '360/Nx', num2str(this.basisCount*180);
                    };
                T = cell2table(SCALES,'VariableNames',{'Scale','Constant','Largest','Smallest'});
                fprintf('\nBasis resolution limits (using Nx patches of size dx):\n');
                disp(T);
                disp(' ');
            end
            % Print information about a particular mode:
            if ~isempty(opts.Results.mode)
                n = opts.Results.mode(:);
                Data = [n this.basisCount./n 2*pi*n/this.basisCount];
                VarNames = {'Mode (-)', 'Wavelength (dx)', 'Wavenumber (1/dx)'};
                fprintf(1, '%s\t\t%s\t\t%s\n', VarNames{:});
                fprintf(1, '%.5g\t\t\t\t%.5g\t\t\t\t\t%.5g\n', Data');
            end
        end
        %% Fourier eigenvalues
        function [eigenvals,wavenumbers,energies,eigamps,dofamps,...
                eigvecmat] = getFourierFootprint(this,varargin)
            % Returns all eigenvalues of the (dimensionless) residual
            % operator in Fourier space for the inviscid advection
            % equation. Look up 'modified wavenumber analysis' for 
            % details; e.g. Van den Abeele (2009) and Alhawwary and Wang
            % (2018).
            %
            % Input (optional)
            %  upwind: amount of upwinding (usually 1 or 0)
            %  wavenumbers: 1D array of (exact, real) wavenumbers associated
            %               with the columns of the eigenvalue matrix
            %               (default: 61 wavemodes, uniformly around zero).
            %  criterion: group eigenmodes such that:
            %                1) each mode's Block wave is smooth
            %             or
            %                2) each mode has similar energy content
            % Output
            %  eigenvals: 2D array of eigenvalues (row: eigenmode; column: 
            %             wavemode). Sorted such that the 1st row is the
            %             "physical" eigenmode.
            %  wavenumbers: actual array of wavenumbers used.
            %
            %  energies: 2D array of relative energies associated to each
            %            eigenmode, for every wavemode.
            %
            %  eigamps: 2D array of eigenmode Fourier coeffs.
            %
            %  dofamps: 2D array of dof-wise Fourier coeffs, as obtained by
            %           projecting the I.C. (== V*eigamps iff t = 0).
            %
            %  eigvecmat: 3D array of eigenvectors, page: wavenumber;
            %               column: eigenvector; row: eigenmode.
            %
            p = inputParser;
            addParameter(p,'upwind',1,@isscalar);
            addParameter(p,'elementCount',500/this.basisCount,@isfinite); % coarse sampling
            addParameter(p,'wavenumbers',[],@isnumeric); % default case depends on 'elementCount' (see below)
            addParameter(p,'criterion','default',@ischar); % how to decide which eigenvalue belongs to which eigenmode
            parse(p,varargin{:});
            % If no wavenumbers are requested, generate from element count:
            if any(strcmp('wavenumbers',p.UsingDefaults))
                M = floor(this.basisCount*p.Results.elementCount/2);
                wavenumbers = 2*pi/p.Results.elementCount*(-M:M);
            else
                wavenumbers = p.Results.wavenumbers;
            end
            % Preallocation:
            eigenvals = complex(zeros(this.basisCount,numel(wavenumbers)));
            energies = zeros(size(eigenvals));
            eigamps = eigenvals;
            dofamps = eigenvals;
            eigvecmat = repmat(permute(eigenvals,[3 1 2]),this.basisCount,1,1);
            % Operator assembly (vectorized):
            [E,leftE,rightE] = this.getFourierMatrices(p.Results.upwind);
            coefsL = exp(-1i*wavenumbers.');
            coefsR = coefsL';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fourier 'interpolatory' coefficients:
            % dofamps = exp(1i*wavenumbers*.5*(1+this.dofCoords));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eigenvalues and eigenvectors:
            for n = 1:numel(wavenumbers)
                R = this.massMatrix \ (E + coefsL(n)*leftE + coefsR(n)*rightE);
                [eigvecmat(:,:,n),D] = eigs(R,this.basisCount);
                eigenvals(:,n) = diag(D);
                % Fourier 'L2' coefficients:
                auxElement = Element(this,[-1 1]); % auxiliary element, used for projecting I.C.
                this.project(auxElement,@(xi) exp(1i*wavenumbers(n)*.5*(1+xi)));
                dofamps(:,n) = auxElement.states;
                % Fourier eigenmode amplitudes:
                eigamps(:,n) = eigvecmat(:,:,n) \ dofamps(:,n);
                % Relative eigenmode energies:
                energies(:,n) = abs(eigamps(:,n)).^2;
                energies(:,n) = energies(:,n)./sum(energies(:,n));
                % Enforce a consistent ordering:
                switch p.Results.criterion
                    case 'none'
                        % do nothing
                    case 'energy' % place the most energetic eigenmode first (screws up the rest)
                        [energies(:,n),i] = sort(energies(:,n),'descend');
                        eigenvals(:,n) = eigenvals(i,n);
                        eigamps(:,n) = eigamps(i,n);
                        eigvecmat(:,:,n) = eigvecmat(:,i,n);
                    otherwise % seek eigenmodes that are smooth, regardless of energy content
                        if n > 1
                            [~,i] = min(abs(eigenvals(:,n) - eigenvals(:,n-1).'),[],2); % i(l): eigenmode that the l-th eigenvalue belongs to
                            j = find(sum(abs(eigenvals(i,n-1) - eigenvals(i,n-1).') < 1e-12) > 1); % ambiguous eigenvalues (assigned to the same eigenmode or almost equal in the previous wavemode)
                            if ~isempty(j) && n > 2
                                k = unique([i(j)' find(sum(1:numel(i) == i) == 0)]); % suspicious eigenmodes (assigned to >1 or <1 eigenvalues, or with almost equal eigenvalues in the previous wavemode)
                                [~,ii] = min(abs(eigenvals(j,n).' - 2*eigenvals(k,n-1) + eigenvals(k,n-2)),[],2); % suspicious eigenmode that each ambiguous eigenvalue belongs to
                                ii(sum(ii == ii') > 1) = find(sum(1:numel(ii) == ii) ~= 1); % assign remaining eigenmodes to remaining eigenvalues 'as is' (if any are still ambiguous)
                                i(j) = k(ii); % j-th eigenvalue actually corresponds to the k(ii)-th eigenmode
                            end
                            eigenvals(i,n) = eigenvals(:,n);
                            energies(i,n) = energies(:,n);
                            eigamps(i,n) = eigamps(:,n);
                            eigvecmat(:,i,n) = eigvecmat(:,:,n);
                        end
                end
            end
            % Sort eigenmodes:
            if n == 1 || strcmp(p.Results.criterion,'none')
                return
            end
            isInRange = wavenumbers/this.basisCount > -pi/2 & wavenumbers/this.basisCount < pi/2; % subset of wavenumbers to consider
            [~,ids] = sort(vecnorm(imag(eigenvals(:,isInRange)) + wavenumbers(:,isInRange),2,2)); % increasing L2-error with exact dispersion relation
            eigenvals = eigenvals(ids,:); % physical first
            energies = energies(ids,:);
            eigamps = eigamps(ids,:);
            eigvecmat = eigvecmat(:,ids,:);
            [~,ids] = sort(polyarea(real(eigenvals(2:end,:)),imag(eigenvals(2:end,:)),2)); % increasing "shadow size" in the complex plane
            eigenvals(2:end,:) = eigenvals(ids+1,:); % "weirdest ones" next
            energies(2:end,:) = energies(ids+1,:);
            eigamps(2:end,:) = eigamps(ids+1,:);
            eigvecmat(:,2:end,:) = eigvecmat(:,ids+1,:);
        end
        %% Modified wavenumbers
        function [k0,varargout] = getModifiedWavenumbers(this,varargin)
            % Returns the exact wavenumbers, then modified ones (one mode
            % per output) of this basis; dimensionless, but not scaled
            % (range is -pi*J to pi*J).
            %
            if nargout > this.basisCount+1
                error('Too many eigenmodes were requested.')
            end
            % Compute the stuff:
            [z,k0] = this.getFourierFootprint(varargin{:});
            % Determine which modes to output:
            varargout(1:nargout-1) = num2cell(1i*z(1:nargout-1,:),2);
        end
        %% Energy-weighted mean modified wavenumber
        function [k0,k] = getMeanWavenumbers(this,varargin)
            % Returns the exact and energy-weighted mean wavenumbers of
            % this basis; dimensionless, but not scaled (range is 0 to
            % pi*J).
            %
            if nargout > this.basisCount+1
                error('Too many eigenmodes were requested.')
            end
            % Compute the stuff:
            [k,k0,e] = this.getFourierFootprint(varargin{:},'criterion','energy');
            k = sum(1i*k.*e,1); % weighted average
        end
        %% Combined-mode modified wavenumber
        function [k0,k,t] = getCombinedWavenumbers(this,varargin)
            % Returns the exact and combined-mode wavenumbers of
            % this basis; dimensionless, but not scaled (range is 0 to
            % pi*J).
            %
            % Arguments
            %  t: (optional) instant(s) when to evaluate the combined mode
            %     phase lag angles and amplification factors.
            %  <...>: passed on to 'Basis.getAngAmp'.
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'t',1,@isnumeric); % t = 1 makes the most sense
            addParameter(p,'mode','combined',@ischar);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            t = p.Results.t(:);
            % Compute phase shift and amplif. factors:
            [k0,angs,amps] = this.getAngAmp(t,varargin{:},'mode','combined');
            % Infer an "effective modified wavenumber":
            k = -angs./t + k0 + 1i*(real(log(amps))./t);
        end
        %% Superconvergence
        function orders = getOrder(this,varargin)
            % Returns a theoretical global order of convergence
            % in a spectral sense, i.e. associated with the dissipation &
            % dispersion error (and not the actual truncation error), a la
            % Vincent et. al. 2011, eq. 4.17.
            %
            % Optional arguments
            %  'k': wavenumbers (not scaled) to be sampled. The output will
            %       be a 2D array, where each row represents an eigenmode
            %       and each column a (given) wavenumber. Eigenmodes are
            %       NOT guaranteed to be sorted.
            %
            %  'upw': upwind ratio for the numerical flux.
            %
            %  'allWavemodes': if true, output every wavemode's order.
            %
            %  'allEigenmodes': if true, output every eigenmode's order.
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'wavenumbers',...
                linspace(-3,3,64)*this.basisCount,... avoids 0 and +-pi*J
                @isnumeric);
            addParameter(p,'allWavemodes',false,@islogical);
            addParameter(p,'allEigenmodes',false,@islogical);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            k0 = p.Results.wavenumbers/2; % coarse mesh wavenumbers
            [k,~,ik] = unique([p.Results.wavenumbers k0]); % k(ik) == [p.Results.wavenumbers k0]
            w = 1i*this.getFourierFootprint('wavenumbers',k,varargin{:}); % modified wavenumbers
            errors = abs(w(:,ik) - k(ik)); % rows <-> eigenmodes
            orders = log(errors(:,1:(end/2))) - log(errors(:,(end/2+1):end));
            orders = orders/log(2) - 1;
            if ~p.Results.allEigenmodes
                orders(2:end,:) = [];
            end
            if ~p.Results.allWavemodes && ~isscalar(k0)
                isOut = (k0 < 0) | [isoutlier(diff(orders(1,:))) false];
                orders = max(orders(:,~isOut),[],2);
            end
        end
        %% Well-resolved wavenumber
        function kf = getResolvingWavenumber(this,varargin)
            % Highest (dimensionless, but not scaled) well-resolved
            % wavenumber, based on relative error threshold (Lele, 1992).
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'rtol',1e-2,@isfinite);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            % Function to solve:
            function f = fun(k0)
                [~,k1] = this.getModifiedWavenumbers(varargin{:},...
                    'criterion','energy','wavenumbers',k0);
                f = abs(k1 - k0) - p.Results.rtol*k0;
            end
            % Iterative search:
            k0 = sinspace(0,pi*this.basisCount,-32);
            k0(~islocalmin(abs(fun(k0)))) = []; % nontrivial root seeds
            kf = fzero(@fun,k0(1)); % converge on first nontrivial root
        end
        %% Cutoff wavenumber
        function kc = getCutoffWavenumber(this,varargin)
            % Estimates a cutoff wavenumber (dimensionless but not scaled)
            % of the physical eigenmode of this basis, using the '1% rule'
            % (Moura et al., 2015).
            %
            % Function to solve:
            function f = fun(k0)
                [~,kM] = this.getModifiedWavenumbers(varargin{:},...
                    'criterion','energy','wavenumbers',k0);
                f = exp(imag(kM)/this.basisCount) - .99;
            end
            % Iterative search:
            kc = abs(fzero(@fun,0)); % root closest to zero (either side)
        end
        %% Dispersion to dissipation ratio
        function [r,k0,k1] = getDispDissRatios(this,mode,varargin)
            % Returns the ratio between dispersion and dissipation effects,
            % a la Adams et al., 2015 (i.e. using group velocity).
            % Only for nonnegative wavenumbers.
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'eps',1e-3,@isfinite);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            switch mode
                case 'single'
                    [k0,k1] = this.getModifiedWavenumbers(varargin{:});
                case 'combined'
                    [k0,k1] = this.getCombinedWavenumbers(varargin{:});
                case 'mean'
                    [z,k0,e] = this.getFourierFootprint(varargin{:});
                    k1 = 1i*z(:,k0 >= 0);
                    e(:,k0 < 0) = [];
                    k0(k0 < 0) = [];
                    r = 0.*k0;
                    for i = 1:this.basisCount
                        num = abs(gradient(real(k1(i,:)),k0) - 1) + p.Results.eps;
                        den = -imag(k1(i,:))/this.basisCount + p.Results.eps;
                        r = r + e(i,:).*num./den;
                    end
                    k1 = sum(k1.*e,1);
                    return
                otherwise
                    error("Unknown 'mode' value '%s'; options are: 'single', 'mean' or 'combined'.",p.Results.mode)
            end
            k1(k0 < 0) = [];
            k0(k0 < 0) = [];
            r = -imag(k1)/this.basisCount + p.Results.eps;
            r = (abs(gradient(real(k1),k0)-1) + p.Results.eps)./r;
        end
        %% Angle and amplification
        function [k,angs,amps,mode] = getAngAmp(this,varargin)
            % Returns the "standard" dissipation and dispersion error
            % measures a la Vanharen et al. (2017): rate of energy loss
            % (amplification factors) and phase lag (angles),
            % respectively, for each wavemode.
            %
            % Arguments
            %  t: (required) array of sample instants (I.C. is at t = 0).
            %  mode: (optional) which eigenmodes to compute the errors for;
            %        valid values are: 'combined' (default), 'mean' (not
            %        recommended), 'physical', 'dominant' and 
            %        'superposition'.
            %
            %  <...>: passed on to Basis.getFourierFootprint.
            %
            % Outputs
            %  k: 1D array of wavenumbers
            %  angs: 2D array of phase lag angles (row: time instant;
            %        column: wavemode).
            %  amps: idem, for amplification factors.
            %  mode: value of the optional input used in the calculation.
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'t',@isnumeric);
            addParameter(p,'mode','combined',@ischar);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            t = p.Results.t(:); % ensure it is a row array
            mode = p.Results.mode; % will output later
            switch mode
                case 'combined'
                    [z,k,~,dofs,dofsL2,V] = this.getFourierFootprint(varargin{:});
                    % Evaluate Fourier coeffs. at given times:
                    t = permute(t,[3 2 1]);
                    dofsL2 = dofsL2.*exp(-1i*k.*t); % dofs x wavemodes x time
                    dofs = dofs.*exp(z.*t); % eigenmodes x wavemodes x time
                    % Convert eigenmode coeffs. to dof coeffs.:
                    dofs = permute(dofs,[1 3 2]); % eigenmodes x time x wavemodes
                    for i = 1:size(V,3)
                        dofs(:,:,i) = V(:,:,i)*dofs(:,:,i); % dofs x time x wavemodes
                    end
                    % Evaluate errors:
                    [angs,amps,ampsL2] = this.getComplexScalarProduct(...
                        [1 2; 1 1; 2 2],... <dofs,dofsL2>, <dofs,dofs>, <dofsL2,dofsL2>
                        permute(dofs,[2 3 1]),... time x wavemodes x dofs
                        permute(dofsL2,[3 2 1])... idem
                    );
                    % Phase angles ('unwrapped' and 'offset-corrected'):
                    angs = angle(angs);
                    [~,id] = min(abs(k));
                    offset = angs(:,id);
                    angs = unwrap(angs,pi,2);
                    offset = offset - angs(:,id);
                    angs = angs + offset;
                    % Amplification factors:
                    amps = sqrt(real(amps./ampsL2)); % remove possible imaginary noise
                    return % finished
                case 'superposition'
                    % Pseudo-combined effect via weighted addition of
                    % individual eigenmodes (Asthana & Jameson's idea)
                    [z,k,e] = this.getFourierFootprint(varargin{:});
                    angs = zeros(numel(t),numel(k));
                    amps = angs;
                    for i = 1:this.basisCount
                        angs = angs + e(i,:).*(imag(z(i,:)) + k).*t;
                        amps = amps + e(i,:).*exp((real(z(i,:)) - 0*k).*t);
                    end
                    return % finished
                case 'mean'
                    [k,kM] = this.getMeanWavenumbers(varargin{:});
                    warning("'Mean' eigenmode tends to overestimate dissipation. Perhaps you meant 'superposition'?")
                case 'physical'
                    [k,kM] = this.getModifiedWavenumbers(varargin{:},'criterion','default');
                case 'dominant'
                    [k,kM] = this.getModifiedWavenumbers(varargin{:},'criterion','energy');
                otherwise
                    error("Invalid value '%s' for argument 'eigenmodes'. Options are: 'combined', 'mean', 'physical', 'dominant' and 'superposition'.",mode)
            end
            % Finish single mode cases:
            angs = (-real(kM) + k).*t;
            amps = exp((imag(kM) - 0*k).*t);
        end
        %% Ang. to amp. ratio norm
        function R = getAngAmpRatioNorm(this,varargin)
            % Returns the L1 norm of the ratio between dispersion and
            % dissipation, in terms of combined mode errors.
            %
            % Only over the badly-resolved range!
            %
            % Parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'t',@isnumeric);
            addParameter(p,'gradient',true,@islogical);
            addParameter(p,'rtol',1e-2,@isfinite);
            addParameter(p,'AbsTol',1e-6,@isfinite);
            addParameter(p,'RelTol',1e-4,@isfinite);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            % Phase lag seeds:
            kf = this.getResolvingWavenumber('rtol',p.Results.rtol);
            [k0,lag0,~] = this.getAngAmp(p.Results.t,varargin{:});
            lag0(:,k0 < kf) = [];
            k0(:,k0 < kf) = [];
            % Function to integrate:
            function r = fun(k)
                % Sample:
                [~,lag,amp] = this.getAngAmp(p.Results.t,varargin{:},'wavenumbers',k);
                % Unwrap the phase lag (multiply-valued):
                [k1,ids] = unique([k0 k]);
                lag1 = [lag0 lag]; % resample
                lag1 = lag1(:,ids); % sort
                lag1 = unwrap(lag1,pi,2); % unwrap
                % Compute the ratio:
                if p.Results.gradient
                    for i = 1:numel(p.Results.t)
                        lag1(i,:) = gradient(lag1(i,:),k1);
                    end
                else
                    lag1 = lag1./k1;
                end
                r = lag1*(k1' == k); % unsort + downsample (linear trans.)
                r = -this.basisCount*abs(r)./log(amp);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 plot(k,lag1*(k1' == k),'+-',k,amp,'x-',k,r,'*-')  %
                 drawnow limitrate                                 %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            % Adaptive quadrature:
            R = integral(@fun,kf,pi*this.basisCount,...
                'AbsTol',p.Results.AbsTol,...
                'RelTol',p.Results.RelTol,...
                'ArrayValued',numel(p.Results.t) > 1);
            R = R/(pi*this.basisCount - kf);
        end
        %% Block wave 3D plot
        function displayModifiedWavenumbers(this,varargin)
            % Plots the modified wavenumbers of the discretization (a
            % priori approach) as 3D waves.
            [z,k] = this.getFourierFootprint(varargin{:});
            z = z/this.basisCount;
            k = k/this.basisCount;
            plot3(k,-imag(z),real(z),k,k,0*k,'--k')
            xlabel('$\kappa^*/J$','Interpreter','LaTex')
            ylabel('$\Re(\tilde{\kappa}^*)/J$','Interpreter','LaTex')
            zlabel('$\Im(\tilde{\kappa}^*)/J$','Interpreter','LaTex')
            view(110,15)
            axis equal
        end
        %% Single-mode dispersion & dissipation
        function displayDispDiss(this,varargin)
            % Plots the real and complex parts of the physical or dominant
            % eigenmode (depending on optional input).
            %
            % Compute:
            [k0,k1] = this.getModifiedWavenumbers(varargin{:});
            k0 = k0/this.basisCount;
            k1 = k1/this.basisCount;
            % Setup plots:
            subplot(2,1,1)
            if isempty(get(gca,'Children'))
                plot(k0([1 end]),k0([1 end]),'--k','DisplayName','Exact disp.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Re}{J}$','Interpreter','LaTex')
            end
            hold on
            subplot(2,1,2)
            if isempty(get(gca,'Children'))
                plot(k0([1 end]),[0 0],'--k','DisplayName','Exact diss.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Im}{J}$','Interpreter','LaTex')
            end
            hold on
            % Plot disp. and diss.:
            subplot(2,1,1)
            plot(k0,real(k1),'DisplayName',this.getName)
            subplot(2,1,2)
            plot(k0,imag(k1),'DisplayName',this.getName)
            % Finishing touches:
            for i = 1:2
                subplot(2,1,i)
                legend('Location','northwest')
                hold off
            end
        end
        %% Energy-weighted mean dispersion & dissipation
        function displayDispDissMean(this,varargin)
            % Plots the real and complex parts of the energy-weighted
            % average among all eigenmodes.
            %
            % IMPORTANT: this is NOT representative of the true spectral
            %            performance of the basis, in general.
            %
            % Compute:
            [k,kMean] = this.getMeanWavenumbers(varargin{:});
            k = k/this.basisCount;
            kMean = kMean/this.basisCount;
            % Setup plots:
            subplot(2,1,1)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),k([1 end]),'--k','DisplayName','Exact disp.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Re}{J}$','Interpreter','LaTex')
            end
            hold on
            subplot(2,1,2)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),[0 0],'--k','DisplayName','Exact diss.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Im}{J}$','Interpreter','LaTex')
            end
            hold on
            % Plot disp. and diss.:
            subplot(2,1,1)
            plot(k,real(kMean),'DisplayName',[this.getName ' (averaged)'])
            subplot(2,1,2)
            plot(k,imag(kMean),'DisplayName',[this.getName ' (averaged)'])
            % Finishing touches:
            for i = 1:2
                subplot(2,1,i)
                legend('Location','northwest')
                hold off
            end
        end
        %% Combined modified wavenumber
        function displayDispDissCombined(this,varargin)
            % Plots the real and imaginary parts of an "effective" modified
            % wavenumber obtained from a combined mode analysis of all
            % eigenmodes.
            %
            % Arguments
            %  <...>: passed on to 'Basis.getCombinedWavenumbers'.
            %
            % Get the stuff:
            [k,kM,t] = this.getCombinedWavenumbers(varargin{:});
            % Proceed as usual:
            k = k/this.basisCount;
            kM = kM/this.basisCount;
            % Setup plots:
            subplot(2,1,1)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),k([1 end]),'--k','DisplayName','Exact disp.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Re}{J}$','Interpreter','LaTex')
            end
            hold on
            subplot(2,1,2)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),[0 0],'--k','DisplayName','Exact diss.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Im}{J}$','Interpreter','LaTex')
            end
            hold on
            % Plot disp. and diss.:
            subplot(2,1,1)
            h = plot(k,real(kM));
            set(h,{'DisplayName'},compose('%s (%s, t* = %g)',this.getName,'combined',t(:)))
            subplot(2,1,2)
            h = plot(k,imag(kM));
            set(h,{'DisplayName'},compose('%s (%s, t* = %g)',this.getName,'combined',t(:)))
            % Finishing touches:
            for i = 1:2
                subplot(2,1,i)
                legend('Location','northwest')
                hold off
            end
        end
        %% Dispersion & dissipation colored by energy content
        function displayDispDissEnergy(this,varargin)
            % Plots the real and complex parts of every modified wavenumber
            % associated to each baseline one, each in a separate scatter
            % subplot.
            %
            % Colors each point according to its relative energy content.
            %
            % Compute:
            [z,k,e] = this.getFourierFootprint(varargin{:});
            z = z/this.basisCount;
            k = k/this.basisCount;
            % Setup plots:
            sz = 1+10*(1-tanh(this.basisCount/10)); % ensures 10 >= sz > 1
            colormap jet
            subplot(2,1,1)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),k([1 end]),'--k','DisplayName','Exact disp.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Re}{J}$','Interpreter','LaTex')
            end
            hold on
            subplot(2,1,2)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),[0 0],'--k','DisplayName','Exact diss.')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{\tilde{\kappa}^*_\Im}{J}$','Interpreter','LaTex')
            end
            hold on
            % Plot disp. and diss.:
            for i = 1:this.basisCount
                subplot(2,1,1)
                scatter(k,-imag(z(i,:)),sz,e(i,:),'filled','DisplayName',sprintf('%s (mode %d)',this.getName,i))
                subplot(2,1,2)
                scatter(k,real(z(i,:)),sz,e(i,:),'filled','DisplayName',sprintf('%s (mode %d)',this.getName,i))
            end
            % Finishing touches:
            for i = 1:2
                subplot(2,1,i)
                legend('Location','northwest')
                hcb = colorbar('eastoutside');
                hcb.Label.String = 'Relative energy';
                hold off
            end
        end
        %% Dispersion, dissipation and their ratio
        function displayDispDissRatio(this,varargin)
            % Plots dispersion, dissipation and their ratio, each in a
            % different subplot, along with the well-resolved range and
            % cutoff wavenumber, in the current figure.
            %
            % Parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'showCutoff','both',@ischar);
            addParameter(p,'mode','combined',@ischar);
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            % Get data:
            [r,k0,k1] = this.getDispDissRatios(p.Results.mode,varargin{:});
            if any(strcmpi(p.Results.showCutoff,{'resolving','DNS','both'}))
                kf = this.getResolvingWavenumber(varargin{:});
            else
                kf = nan;
            end
            if any(strcmpi(p.Results.showCutoff,{'1%','LES','both'}))
                kc = this.getCutoffWavenumber(varargin{:});
            else
                kc = nan;
            end
            % Scale it:
            k0 = k0/this.basisCount;
            k1 = k1/this.basisCount;
            kf = kf/this.basisCount;
            kc = kc/this.basisCount;
            % Set it up:
            refData = {k0([1 end]),[0 0],[1 1]};
            numData = {real(k1) imag(k1) r};
            yLabels = {
                '$$\frac{\tilde{\kappa}^*_\Re}{J}$$'
                '$$\frac{\tilde{\kappa}^*_\Im}{J}$$'
                '$$\frac{J |\frac{d\tilde{\kappa}^*_\Re}{d\kappa} - 1|}{-\tilde{\kappa}^*_\Im}$$'
            };
            % Plot it: 
            for i = 1:3
                subplot(3,3,3*(i-1)+[1 2])
                if isempty(get(gca,'Children'))
                    plot(k0([1 end]),refData{i},'--k','DisplayName','Ideal')
                end
                hold on
                h = plot(k0,numData{i},'DisplayName',this.getName);
                plot([kf kf],[min(numData{i}) max(numData{i})],'--','Color',h.Color,'DisplayName',['Resolv. cutoff; ' this.getName])
                plot([kc kc],[min(numData{i}) max(numData{i})],':','Color',h.Color,'DisplayName',['1% diss. cutoff; ' this.getName])
                hold off
                xlabel('$$\frac{\kappa^*}{J}$$','Interpreter','Latex')
                ylabel(yLabels{i},'Interpreter','Latex')
            end
            % Add a separate legend:
            g = subplot(3,3,3:3:9);
            set(g,'Visible','off')
            h = get(gcf,'Children');
            legend(g,h(end).Children,'Location','best')
        end
        %% Alternative disp./diss. ratio
        function displayAngAmpRatio(this,t,varargin)
            % Plots disp./diss. ratios resulting from angle shift and
            % amplification factor at requested time instants.
            %
            % Arguments:
            %  t: (required) time(s) at which to sample the ratio.
            %  eps: (optional) custom small value to avoid indetermination.
            %  <...>: passed on to 'basis.getAngAmp'.
            %
            % Parse:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'eps',1e-3,@isfinite);
            addParameter(p,'showCutoff','neither',@ischar); % unused
            parse(p,varargin{:});
            varargin = [fields(p.Unmatched) struct2cell(p.Unmatched)]';
            % Get data:
            [k,angs,amps,mode] = this.getAngAmp(t,varargin{:});
            for i = numel(t):-1:1
                r(i,:) = gradient(angs(i,:),k);
            end
            r = (abs(r) + p.Results.eps)./(-log(amps)/this.basisCount + p.Results.eps);
            % Set it up:
            k = k/this.basisCount;
            refData = {[0 0],[1 1],[1 1]};
            numData = {angs/this.basisCount amps r};
            yLabels = {
                '$$\frac{\Delta\Psi^*}{J}$$'
                '$$G$$'
                '$$\frac{|\frac{d(\Delta\Psi^*)}{d\kappa^*}|}{\frac{-\ln G}{J}}$$'
            };
            % Plot it: 
            for i = 1:3
                subplot(3,3,3*(i-1)+[1 2])
                if isempty(get(gca,'Children'))
                    plot(k([1 end]),refData{i},'--k','DisplayName','Ideal')
                end
                hold on
                h = plot(k,numData{i});
                set(h,{'DisplayName'},compose('%s (%s, t = %g)',this.getName,mode,t(:)));
                hold off
                xlabel('$$\frac{\kappa^*}{J}$$','Interpreter','Latex')
                ylabel(yLabels{i},'Interpreter','Latex')
            end
            % Add a separate legend:
            g = subplot(3,3,3:3:9);
            set(g,'Visible','off')
            h = get(gcf,'Children');
            legend(g,h(end).Children,'Location','best')
        end
        %% Eigenmode energy content
        function displayEnergy(this,varargin)
            % Plots the relative energy content in each eigenmode.
            %
            [~,k,e] = this.getFourierFootprint(varargin{:});
            hold on
            h = plot(k/this.basisCount,e);
            set(h,{'DisplayName'},...
                compose('%s (mode %d)',this.getName,1:this.basisCount)')
            hold off
            xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
            ylabel('$\Gamma$','Interpreter','LaTex')
            legend('Location','west')
        end
        %% Combined mode plots
        function displayAngAmp(this,t,varargin)
            % Plots the disp. and diss. errors as in a combined mode
            % analysis a la Alhawwary and Wang (2018).
            %
            % Arguments
            %  t: time instant(s) for which to compute the errors
            %  <...>: passed on to Basis.getDispDissErrors
            %
            % Get the data:
            [k,angs,amps,mode] = this.getAngAmp(t,varargin{:});
            k = k/this.basisCount;
            % Setup plots:
            subplot(2,1,1)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),[0 0],'--k','DisplayName','Exact (no phase shift)')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$\frac{|\Delta\psi(t^*)|}{J}$','Interpreter','LaTex')
            end
            hold on
            subplot(2,1,2)
            if isempty(get(gca,'Children'))
                plot(k([1 end]),[1 1],'--k','DisplayName','Exact (no damping)')
                xlabel('$\frac{\kappa^*}{J}$','Interpreter','LaTex')
                ylabel('$G(t^*)$','Interpreter','LaTex')
            end
            hold on
            % Plot disp. and diss.:
            subplot(2,1,1)
            h = plot(k,abs(angs)/this.basisCount);
            set(h,{'DisplayName'},compose('%s (%s, t* = %g)',this.getName,mode,t)')
            set(gca,'YScale','log')
            subplot(2,1,2)
            h = plot(k,amps);
            set(h,{'DisplayName'},compose('%s (%s, t* = %g)',this.getName,mode,t)')
            % Finishing touches:
            for i = 1:2
                subplot(2,1,i)
                legend('Location','best')
                xlim([0 pi])
                hold off
            end
        end
        %% Combined wavenumbers vs. t
        function displayCombinedWavenumberSensitivity(this,varargin)
            % This function plots the vector 2-norm of the combined-mode
            % wavenumber as a function of time.
            %
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'t',logspace(-10,2,100),@isnumeric);
            parse(p,varargin{:});
            [k,kM,t] = this.getCombinedWavenumbers(varargin{:},'t',p.Results.t);
            aux = vecnorm(kM,2,2);
            aux = gradient(aux,k);
            hold on
            plot(t,aux,'+-','DisplayName',this.getName)
            hold off
            set(gca,'XScale','log')
            xlabel('$$t^*$$','Interpreter','latex')
            ylabel('$$\frac{d \vert\vert \tilde{\kappa}^* \vert\vert}{d\kappa^*}$$','Interpreter','latex')
            legend('Location','Best')
        end
    end
    methods (Access = protected)
        %% Assemble Fourier residual matrices (DG)
        function [A,B,C] = getFourierMatrices(this,beta)
            % Base implementation (valid for all pure DG methods; not FR).
            A = 2*this.gradientMatrix.' + sparse((1-beta)*this.left*this.left.' - (1+beta)*this.right*this.right.');
            B = (1+beta)*sparse(this.left*this.right.'); % upwind
            C = (-1+beta)*sparse(this.right*this.left.'); % downwind
        end
        %% Complex scalar product(s)
        function [varargout] = getComplexScalarProduct(this,perms,varargin)
            % Computes the scalar product operator between two complex
            % functions assumed to exist in the trial space spanned by
            % this basis.
            %
            % This is the operator that defines the L2 norm in a complex
            % space (or something like that, I'm not a mathematician); in
            % any case, it is the one defined in Vanharen et al. (2017) and
            % used in Alhawwary and Wang (2018) to study combined mode
            % spectral analysis stuffs.
            %
            % This version assumes that this basis spans both trial and
            % test spaces; for e.g. FR, it needs to be overriden.
            %
            % Arguments
            %  perms: 2D array of permutations to compute. Row: output
            %         variable (permutation); column: input variable.
            %         Number of columns must be 2 (this is a binary op.),
            %         and the number of rows must be <=(nargin-2)^2.
            %  varargin: list of 3D arrays of dofs; each page is assumed to
            %            correspond to a different basis function. The
            %            other two are reserved for, e.g., time instants
            %            and wavenumbers.
            %
            massVec = permute(... array of nonzero mass matrix entries
                full(this.massMatrix(this.pairs(3,:))),[1 3 2]);
            % Loop over requested permutations:
            for i = nargout:-1:1
                varargout{i} = sum(... discrete complex inner product
                    varargin{perms(i,1)}(:,:,this.pairs(1,:)).*...
                    conj(varargin{perms(i,2)}(:,:,this.pairs(2,:))).*...
                    massVec,3); % unrolled AND vectorized
            end
        end
    end
    methods (Static, Access = protected, Sealed)
        %% Default type (Override)
        function default_object = getDefaultScalarElement
            % When an array of basis children is preallocated (e.g.
            % bases(2) = DGSEM) any non-specified elements are filled
            % with DG instances (e.g. class(bases(1)) == 'DG').
            default_object = DG;
        end
    end
end