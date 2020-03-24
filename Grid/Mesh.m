classdef Mesh < handle
    % Class that groups an arbitrary number of elements/patches and couples
    % them with a (possibly different) number of unique bases. Interface to
    % many mesh-wide methods and properties.
    %
    properties
        maxDegree
        minDegree
        maxBasisCount
        minBasisCount
        dofCount % total number of degrees of freedom, per equation
        edgeCount
        edges
        elementCount
        elements
        % Row array of finite-dimensional spaces (in reference element
        % coordinates) in which the approximate solution lives
        bases
    end
    methods
        %% Constructor
        function this = Mesh(varargin)
            % Intantiates a mesh instance according to parsed name-value
            % argument pairs. All arguments are self-explanatory.
            %
            % Check:
            if nargin < 2
                return % default constructor
            end
            % Input parser:
            p = inputParser;
            addRequired(p,'bases',@this.validate_bases)
            addRequired(p,'edgeCoords',@this.validate_edgeCoords)
            addOptional(p,'elementCount',numel(varargin{2})-1,@this.validate_elementCount)
            addOptional(p,'edgeDistribution',"uniform",@this.validate_edgeDistribution)
            addOptional(p,'clusteringFactor',0,@this.validate_clusteringFactor)
            addParameter(p,'degrees',[varargin{1}.degree],@this.validate_degrees)
            parse(p,varargin{:});
            % Initialize array of distinct bases:
            if numel(p.Results.bases) == 1 % smallest set of bases with specified degrees
                [degrees,~,element2basis] = unique(p.Results.degrees);
                this.bases = arrayfun(@(x) p.Results.bases.clone(x),degrees);
            else % specified set of bases (possibly repeated)
                element2basis = 1:numel(p.Results.bases);
                if numel(p.Results.degrees) == 1 % enforce specified degree for all bases
                    this.bases = arrayfun(@(x) x.clone(p.Results.degrees),p.Results.bases);
                elseif numel(p.Results.bases) == numel(p.Results.degrees) % overwrite each basis' degree
                    this.bases = p.Results.bases;
                    for k = 1:numel(this.bases)
                        this.bases(k) = p.Results.bases(k).clone(p.Results.degrees(k));
                    end
                else % ignore degrees
                    warning('Specified degrees will be ignored.')
                    this.bases = p.Results.bases;
                end
            end
            % Initialize array of edge coordinates:
            this.edgeCount = numel(p.Results.edgeCoords);
            this.elementCount = p.Results.elementCount;
            if this.edgeCount ~= this.elementCount+1 % generate edge coordinates
                this.edgeCount = p.Results.elementCount + 1;
                switch this.validate_edgeDistribution(p.Results.edgeDistribution)
                    case "uniform"
                        x = linspace(p.Results.edgeCoords(1),p.Results.edgeCoords(end),this.edgeCount);
                    case "cosine"
                        x = cosspace(p.Results.edgeCoords(1),p.Results.edgeCoords(end),this.edgeCount,p.Results.clusteringFactor);
                    case "sine"
                        x = sinspace(p.Results.edgeCoords(1),p.Results.edgeCoords(end),this.edgeCount,p.Results.clusteringFactor);
                    case "logarithmic"
                        x = logspace(p.Results.edgeCoords(1),p.Results.edgeCoords(end),this.edgeCount);
                end
            else % use given coordinates (ignore distribution)
                x = p.Results.edgeCoords;
            end
            % Initialize array of elements (with circular repetition, if necessary):
            this.elements = arrayfun(@Element,x(1:end-1),x(2:end),this.bases(element2basis(mod(0:this.elementCount-1,numel(element2basis))+1)));
            % Initialize cell array of element edges:
            this.edges = {LeftBoundary(x(1),this.elements(1)) arrayfun(@Edge,x(2:end-1),this.elements(1:end-1),this.elements(2:end)) RightBoundary(x(end),this.elements(end))};
            % Initialize all remaining properties:    
            this.minDegree = min([this.bases.degree]);
            this.maxDegree = max([this.bases.degree]);
            this.minBasisCount = min([this.bases.basisCount]);
            this.maxBasisCount = max([this.bases.basisCount]);    
            this.dofCount = sum([this.elements.dofCount]);
        end
        %% Discretization summary
        function info = getInfo(this)
            % Provides concise information about the discretization of the
            % domain and the solution space (in one line).
            info = this.bases(1).getName;
            if numel(this.bases) > 1
                info = sprintf('%s... (%d more), p \\in [%d,%d]',info,numel(this.bases)-1,this.minDegree,this.maxDegree);
            end
            info = sprintf('%s; N_\\Omega = %d, N = %d',info,this.elementCount,this.dofCount);
        end
        %% Compute mesh residuals
        function computeResiduals(this,physics)
            % Updates the residuals of all cells in a mesh, using the
            % spatial discretization scheme assigned to the mesh. 
            %
            % Four-stage process:
            % 1) Computes the state at both edges of each element.
            % 2) Computes the Riemann flux at each interior interface; 
            %    assigns it (conventional sign) to each adjacent element.
            % 3) Computes Riemann flux at domain boundaries (using boundary
            %    conditions).
            % 4) Applies the spatial discretizationn operator to each
            %    element.
            %
            for element = this.elements
                element.localTimeDelta = inf;
                element.interpolateStateAtEdges;
            end
            for edge = this.edges{2:end-1} % interior edges only
                edge.computeFlux(physics);
            end
            physics.applyBoundaryConditions(this); % boundary edges
            for element = this.elements
                element.computeResiduals(physics);
            end
        end
        %% Extract matrices of nodal + edge values of the entire mesh
        function q = getNodalCoeffs(this,eqs)
            % eqs: 1D array of system components to extract.
            % q(i,j,k): cell array of 2D arrays of nodal (and edge) state 
            %           values. 
            %     i: system component
            %     j = 1: left edge value
            %     1 < j < Nb+2: nodal values
            %     j = Nb+2: right edge value
            %     k: element
            %    Elements with Nb < max(Nb) have their j > Nb + 2 components
            %    padded with NaN.
            nEqs = length(eqs);
            q = nan(nEqs,this.maxBasisCount+2,this.elementCount); % NaN padding
            k = 0;
            for element = this.elements
                k = k + 1;
                % Nodal values:
                if element.basis.degree > 0
                    q(:,2:element.basis.basisCount+1,k) = element.states(eqs,:);
                end
                % Edge values:
                element.interpolateStateAtEdges;
                q(:,1,k) = element.stateL;
                q(:,element.basis.basisCount+2,k) = element.stateR;
            end
        end
        %% Sample the states globally
        function states = sample(this,x)
            % x: 1D row array of sample locations
            % samples: 2D array of sampled (column) state vectors at each x
            %  i: state vector component
            %  j: location index
            states = zeros(size(this.elements(1).states,1),length(x));
            for element = this.elements
                ids = x >= element.xL & x <= element.xR;
                y = x(ids); % sample locations within current element
                if ~isempty(y)
                    states(:,ids) = element.interpolateStateAtCoords(y);
                end
            end
        end
        %% Sample the residuals globally
        function residuals = sampleResidual(this,x)
            % x: 1D row array of sample locations
            % samples: 2D array of sampled (column) residual vectors at x
            %  i: state vector component
            %  j: location index
            residuals = zeros(size(this.elements(1).states,1),length(x));
            for element = this.elements
                ids = x >= element.xL & x <= element.xR;
                y = x(ids); % sample locations within current element
                if ~isempty(y)
                    residuals(:,ids) = element.interpolateResidualAtCoords(y);
                end
            end
        end
        %% Sample DOF quantities globally
        function varargout = sampleDOFs(this,x,varargin)
            % Retrieve any quantity defined at the DOFs of each element in
            % the mesh. Selection via propery names. No checks are made.
            %
            % Arguments
            %  x: 1D row array of sample locations
            %  varargin: list of property names to sample (e.g. states, residuals, fluxes)
            %  varargout: list of 2D arrays of sampled (column) vectors
            %   i: vector component
            %   j: location index
            N = length(varargin);
            if N == 0
                varargout = {};
                return
            end
            varargout = zeros(size(this.elements(1).states,1),length(x),N);
            for element = this.elements
                ids = x >= element.xL & x <= element.xR;
                y = x(ids); % sample locations within current element
                if ~isempty(y)
                    z = element.mapToReference(y);
                    phi = element.basis.sampleAt(z);
                    for i = 1:N
                        varargout(:,ids,i) = element.(varargin{i})*phi;
                    end
                end
            end
            varargout = num2cell(varargout,[1 2]);
        end
        %% Retrieve element centroid locations globally
        function x = getElementLocations(this)
            % Preallocate, then fill:
            x = zeros(1,this.elementCount);
            n = 0;
            for element = this.elements
                n = n + 1;
                x(n) = 0.5*(element.xL + element.xR);
            end
        end
        %% Retrieve breakpoint locations globally
        function x = getBreakLocations(this,makeUnique)
            % Remove duplicated element edge breakpoints?
            if nargin > 1 && makeUnique
                l = 0; % yes
            else
                l = 1; % no (default)
            end
            N = 1;
            % Count the total number of breakpoints:
            for element = this.elements
                N = N + length(element.basis.breakCoords) - 1;
            end
            % Preallocate, then fill:
            x = zeros(1,N);
            n = 1;
            for element = this.elements
                m = n + length(element.basis.breakCoords) - 1;
                x(n:m) = element.getBreakCoords;
                n = m + l;
            end
        end
        %% Retrieve nodal locations globally
        function x = getNodeLocations(this)
            N = 0;
            % Count the total number of nodes:
            for element = this.elements
                N = N + element.basis.basisCount;
            end
            % Preallocate, then fill:
            x = zeros(1,N);
            n = 1;
            for element = this.elements
                m = n + element.basis.basisCount - 1;
                x(n:m) = element.getNodeCoords;
                n = m + 1;
            end
        end
        %% Retrieve quadrature locations globally
        function x = getGaussLocations(this)
            x = cell(1,this.elementCount);
            i = 1;
            for element = this.elements
                x{i} = reshape(element.getGaussCoords,1,[]); % ensures that coordinates are a row array
                i = i + 1;
            end
            x = cell2mat(x);
        end
        %% Retrieve DOF locations globally
        function x = getDofLocations(this)
            % Coordinates associated with the degrees of freedom of every
            % discretization, i.e. quadrature locations (Lagrange and 
            % Legendre) or control point locations (B-spline).
            x = cell(1,this.elementCount);
            i = 1;
            for element = this.elements
                x{i} = reshape(element.getDofCoords,1,[]);
                i = i + 1;
            end
            x = cell2mat(x);
        end
        %% Compute mass norms
        function mass = getSolutionMass(this,eqs)
            % Computes the "mass norm" (non-normalized L1, without absolute 
            % values) of solution vector components via Gauss quadrature of
            % optimal order for each element.
            %
            % eqs: column array (1 by default)
            if nargin < 2
                eqs = 1; % assume scalar case
            end
            % Preallocate:
            mass = zeros(length(eqs),1);
            dSigma = 0;
            % Loop over elements:
            for element = this.elements
                if element.basis.isHybrid
                    B = element.basis.vandermonde;
                    modes = element.basis.span2mode;
                elseif element.basis.isNodal
                    B = eye(element.basis.basisCount);
                    modes = 1:element.basis.basisCount;
                else
                    B = element.basis.vandermonde;
                    modes = 1:element.basis.basisCount;
                end
                % Loop over knot spans:
                for l = 1:size(B,3)
                    dSigma = dSigma + 2; % assumes length of 2
                    mass = mass +...
                        element.states(eqs,modes(l,:))...
                        *B(:,:,l)*element.basis.gaussWeights;
                end
            end
            % Convert from knot span units to physical domain units:
            mass = mass/dSigma*(this.edges{end}.coord - this.edges{1}.coord);
        end
        %% Total variation in the means
        function tvm = getTVM(this,eqs)
            % Computes the total variation in the (patch-wide) means of 
            % solution vector components via Gauss quadrature of optimal 
            % order.
            %
            % eqs: column array (1 by default)
            if nargin < 2
                eqs = 1; % assume scalar case
            end
            % Preallocate:
            masses = zeros(length(eqs),this.elementCount); % patch-wise mass norms
            % Compute the mass in each patch:
            k = 0;
            for element = this.elements
                k = k + 1;
                if element.basis.isHybrid
                    B = element.basis.vandermonde;
                    modes = element.basis.span2mode;
                elseif element.basis.isNodal
                    B = eye(element.basis.basisCount);
                    modes = 1:element.basis.basisCount;
                else
                    B = element.basis.vandermonde;
                    modes = 1:element.basis.basisCount;
                end
                % Loop over knot spans:
                for l = 1:size(B,3)
                    masses(:,k) = masses(:,k) +...
                        element.states(eqs,modes(l,:))...
                        *B(:,:,l)*element.basis.gaussWeights;
                end
                % Normalize by the length of the patch:
                masses(:,k) = masses(:,k)/(2*size(B,3));
            end
            % Compute the total variation in the means:
            tvm = sum(abs(masses(:,2:end) - masses(:,1:end-1)),2);
        end
        %% Compute solution norms (using adaptive Gauss-Konrod quadrature)
        function norms = getSolutionNorm(this,p,eqs,x)
            % Computes the p-norm, generalized from the L2 norm (a la Wang
            % et al. 2013), of the approximate solution over the entire
            % mesh.
            %
            %  p: p-norm parameter (2 by default)
            %  eqs: column array (1 by default)
            %  x: row array of expected singular locations (optional)
            %
            if nargin < 2
                p = 2;
            end
            if nargin < 3
                eqs = 1;
            end
            if nargin < 4
                x = [];
            end
            % Determine integration sub-intervals:
            subs = [];
            for element = this.elements
                if isa(element.basis,'Bspline')
                    xi = element.basis.breakCoords;
                    xi = element.mapFromReference(xi);
                else
                    xi = [element.xL element.xR];
                end
                if ~isempty(x)
                    xi = sort([xi x(x < xi(end) & x > xi(1))]);
                end
                subs = cat(2,subs,[xi(1:end-1); xi(2:end)]);
            end
            % Compute the norm:
            if isnumeric(p)
                norms = Algorithms.quadvgk(@(x) abs(this.sample(x)).^p,subs,length(eqs));
                norms = nthroot(norms./(this.edges{end}.coord - this.edges{1}.coord),p);
            elseif strcmp(p,'mass')
                norms = Algorithms.quadvgk(@(x) this.sample(x),subs,length(eqs));
            else
                error('Unknown or unsupported norm.')
            end
        end
        %% Compute error norms (using adaptive Gauss-Konrod quadrature)
        function norms = getErrorNorm(this,FUNs,p,eqs,x)
            % Computes the L2 norm of the approximate solution a la Wang et
            % al. (2013).
            %
            % FUNs: function handle
            % p: p-norm parameter (2 by default)
            % eqs: column array (1 by default)
            % x: row array of expected singular locations (optional)
            if nargin < 3
                p = 2;
            end
            if nargin < 4
                eqs = 1; % assume FUNs is a scalar equation
            end
            if nargin < 5
                x = [];
            end
            % Determine integration sub-intervals:
            subs = [];
            for element = this.elements
                if isa(element.basis,'Bspline')
                    xi = element.basis.breakCoords;
                    xi = element.mapFromReference(xi);
                else
                    xi = [element.xL element.xR];
                end
                if ~isempty(x)
                    xi = sort([xi x(x < xi(end) & x > xi(1))]);
                end
                subs = cat(2,subs,[xi(1:end-1); xi(2:end)]);
            end
            norms = Algorithms.quadvgk(@(x) abs(this.sample(x) - FUNs(x)).^p, subs, length(eqs));
            norms = nthroot(norms./(this.edges{end}.coord - this.edges{1}.coord),p);
        end
        function [R,success] = getTotalVariation(this,maxIters)
            % Compute the total variation of the current approximate
            % solution stored in the mesh. Convergence acceleration using
            % linear Richardson extrapolation.
            %
            if nargin == 1
                maxIters = 30;
            end
            success = 0;
            %x = this.getNodeLocations;
            x = this.getGaussLocations; % global quadrature point locations
            q = this.sample(x); % state vector at every quadrature location (initial samples)
            function tv = TV
                tv = sum(abs(diff(q,1,2)),2);
            end
            A(:,2) = TV;
            R0 = inf;
            % Iterate until reaching tolerance or excessive cost:
            for iter = 1:maxIters
                [q,x] = this.resampleBisection(q,x);
                A = horzcat(A(:,2),TV);
                if all(abs(diff(A,1,2)) < 1e-10)
                    R = A(:,end);
                    return
                end
                R = Algorithms.richardsonExtrapolate(A);
                if all(abs(R - R0) < 1e-10)
                    return
                end
                if length(x) > 1e6 || any(isnan(R)) || any(R > 1e200)
                    %warning('Failed to converge on tolerance.')
                    break
                end
                R0 = R;
            end
            success = 1;
        end
        function [q,x] = resampleBisection(this,q0,x0)
            % Refine (in a nested way) a set of solution samples (from
            % this mesh), by adding to the provided sampled states and 
            % locations new ones obtained via bisection. Only the new 
            % samples are computed.
            %
            x = zeros(1,2*length(x0)-1);
            q = zeros(size(q0,1),length(x));
            ids = logical(mod(1:length(x),2)); % old entries
            x(ids) = x0;
            q(:,ids) = q0;
            ids = ~ids; % new entries
            x(ids) = .5*(x0(2:end) + x0(1:end-1));
            q(:,ids) = this.sample(x(ids));
        end
    end
    methods (Static,Access = protected)
        %% Custom validation functions
        function validate_edgeCoords(x)
            % Throws an error if the input is not a valid set of edge 
            % coordinates, i.e. an array of 2 or more finite real values,
            % in increasing order.
            validateattributes(x,{'numeric'},{'nonempty','finite','increasing'})
            if isscalar(x)
                error('Expected input to be non-scalar.')
            end
        end
        function validate_bases(x)
            % Throws an error if the input is not a valid set of bases.
            validateattributes(x,{'Basis'},{'nonempty','vector'})
        end
        function validate_elementCount(x)
            % Throws an error if the input is not a valid number of
            % elements.
            validateattributes(x,{'numeric'},{'integer','scalar'})
        end
        function varargout = validate_edgeDistribution(x)
            % Throws an error of the input is not a valid kind of edge
            % distribution. If an output is requested, returns the full 
            % name of the matching distribution.
            %
            name = validatestring(x,["uniform","cosine","sine","logarithmic"]);
            if nargout > 0
                varargout{1} = name;
            end
        end
        function validate_clusteringFactor(x)
            % Throws an error if input is not a valid clustering factor.
            validateattributes(x,{'numeric'},{'finite','scalar'})
        end
        function validate_degrees(x)
            % Throws an error if the input is not a valid set of degrees.
            validateattributes(x,{'numeric'},{'integer'})
        end
    end
end