classdef Mesh < handle
    properties
        maxDegree
        minDegree
        maxBasisCount
        minBasisCount
        edgeCount
        edges
        elementCount
        elements
        bases % collection of unique finite-dimensional spaces (in
              % reference element coordinates) of a DG sub-type in
              % which the solution will be approximated
        dofCount % total number of degrees of freedom, per equation
        physics
    end
    methods
        %% Constructor
        function mesh = Mesh(x,p,discretization,physics)
            % Arguments
            %  x: element edge locations
            %  p: element-wise polynomial degree
            %  discretization: spatial discretization instance
            if nargin > 0
                mesh.edgeCount = length(x);
                mesh.elementCount = mesh.edgeCount-1;
                if length(p) ~= mesh.elementCount
                    p = p(1)*ones(1,mesh.elementCount);
                end
                mesh.minDegree = min(p);
                mesh.maxDegree = max(p);
                % Initialize collection of discretization instances:
                [p,~,collectionIds] = unique(p);
                mesh.bases = arrayfun(@discretization.clone,p);
                mesh.minBasisCount = min(cell2mat({mesh.bases.basisCount}));
                mesh.maxBasisCount = max(cell2mat({mesh.bases.basisCount}));
                % Initialize array of elements:
                mesh.elements = arrayfun(@Element,x(1:end-1),x(2:end),...
                    mesh.bases(collectionIds'));
                mesh.dofCount = sum(cell2mat({mesh.elements.dofCount}));
                % Set connectivity:
                mesh.edges = arrayfun(@Edge,x(2:end-1),...
                    mesh.elements(1:end-1),mesh.elements(2:end));
                mesh.edges = {
                    LeftBoundary(x(1),mesh.elements(1))...
                    mesh.edges...
                    RightBoundary(x(end),mesh.elements(end))
                    };
                % Store physics handle:
                mesh.physics = physics;
            end
        end
        %% Compute mesh residuals
        function computeResiduals(this)
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
                edge.computeFlux(this.physics);
            end
            this.physics.applyBoundaryConditions(this); % boundary edges
            for element = this.elements
                element.computeResiduals(this.physics);
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
        %% Retrieve element centroid locations
        function x = getElementLocations(this)
            % Preallocate, then fill:
            x = zeros(1,this.elementCount);
            n = 0;
            for element = this.elements
                n = n + 1;
                x(n) = 0.5*(element.xL + element.xR);
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
                if element.basis.isHybrid
                    x{i} = reshape(element.getControlCoords,1,[]);
                else
                    x{i} = reshape(element.getGaussCoords,1,[]);
                end
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
end