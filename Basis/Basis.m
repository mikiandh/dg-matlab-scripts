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
        function [eigenvals,wavenumbers] = getFourierFootprint(this,beta,wavenumbers)
            % Returns all eigenvalues of the (dimensionless) residual
            % operator in Fourier space for the inviscid advection
            % equation. Look up 'modified wavenumber analysis' for 
            % details (e.g. Van den Abeele, 2009).
            %
            % Input
            %  beta: upwind ratio (usually 1 to 0)
            %  wavenumbers: 1D array of (exact, real) wavenumbers associated
            %               with the columns of the eigenvalue matrix
            %               (default: 61 wavemodes, uniformly around zero).
            % Output
            %  eigenvals: 2D array of eigenvalues (row: eigenmode; column: 
            %             wavemode). Sorted such that the 1st row is the
            %             "physical" eigenmode.
            %  wavenumbers: actual array of wavenumbers used.
            %
            if nargin < 2
                beta = 1; % default to upwind
            end
            if nargin < 3
                wavenumbers = pi*this.basisCount*linspace(-1,1,this.basisCount*60); % 60 generating patterns (resolution)
            end
            % Preallocation:
            eigenvals = complex(zeros(this.basisCount,numel(wavenumbers)));
            % Operator assembly (vectorized):
            [E,leftE,rightE] = this.getFourierMatrices(beta);
            coefsL = exp(-1i*wavenumbers.');
            coefsR = coefsL';
            % Eigenvalues:
            for n = 1:numel(wavenumbers)
                R = this.massMatrix \ (E + coefsL(n)*leftE + coefsR(n)*rightE);
                eigenvals(:,n) = eigs(R,this.basisCount);
                % Enforce a consistent ordering:
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
                end
            end
            % Sort eigenmodes:
            isInRange = wavenumbers/this.basisCount > -pi/2 & wavenumbers/this.basisCount < pi/2; % subset of wavenumbers to consider
            [~,ids] = sort(vecnorm(imag(eigenvals(:,isInRange)) + wavenumbers(:,isInRange),2,2)); % increasing L2-error with exact dispersion relation
            eigenvals = eigenvals(ids,:); % physical first
        end
        function displayModifiedWavenumbers(this,varargin)
            % Plots the modified wavenumbers of the discretization (a
            % priori approach) as 3D waves.
            [z,k] = this.getFourierFootprint(varargin{:});
            z = z/this.basisCount;
            k = k/this.basisCount;
            plot3(k,-imag(z),real(z),k,k,0*k,'--k')
            xlabel('$\kappa/J$','Interpreter','LaTex')
            ylabel('$\Re(\tilde{\kappa})/J$','Interpreter','LaTex')
            zlabel('$\Im(\tilde{\kappa})/J$','Interpreter','LaTex')
            view(110,15)
            axis equal
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