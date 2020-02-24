classdef Basis < matlab.mixin.SetGet
    properties (Abstract)
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
        massMatrix % inner products between every two basis functions
        gradientMatrix % inner products between every basis function (row) and first derivative (column)
        nodeCoords
        gaussCoords
        gaussWeights
        pairs % 2D array of test function/basis function pairs within mutual nonzero support; 1st row: test function, r; 2nd row: basis function, j; 3rd row: linear index; 4th row: transposed linear index
        edges % idem, but excludes pairs of repeated control points, i.e. j == r ("virtual" edges a la Kuzmin et al. 2012)
    end
    methods (Abstract)
        clone(prototype,degree) % required for the prototype instantiation pattern
        computeResiduals(this,element,physics) % evaluate discrete spatial operator, i.e residual (a la Wang et al. 2013)
        sampleAt(this,x) % sample basis components at given locations
        interpolate(mesh,limiter,fun) % sample-based projection
        project(mesh,limiter,fun) % L2 projection (conservative)
        %%% integrateFun(fun,q) % numerically integrate a given function over the reference patch
    end
    methods
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
        function x = getGaussCoords(this)
            % Retrieves quadrature point locations in reference patch 
            % coordinates. Default implementation (override as necessary).
            x = this.gaussCoords;
        end
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
    end
end