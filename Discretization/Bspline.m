classdef Bspline < Basis
    properties
        % Flags:
        isNodal = false; % can be set to 'true' for debugging purposes
        isModal = false;
        isHybrid = true; % combines DG and CG elements
        % Miscellanea:
        smoothness % differentiability class for the entire patch (if applicable)
        controlCoords % abscissae of all control points in a patch (in reference patch coordinates)
        breakCoords % abscissae of all breakpoints in a patch (in reference patch coordinates)
        nonzeroSpanCount % number of non-vansihing spans
        % B-splines:
        knots % open knot vector
        knotDiffs % knot differences
        knotCount % number of knots
        % Operators:
        span2knot % connectivity vector between a non-vanishing knot span index and the position of its left knot in the knot vector
        span2mode % connectivity matrix between non-vanishing knot spans (rows) and basis functions with support on each (columns)
        vandermonde % 3D array of basis functions (rows) sampled at Gauss points (columns) of each nonzero knot span (pages)
        derivatives % idem, for the first derivative of the basis functions
        invVandermonde % 3D array of Gauss quadrature point contributions (rows) to modal components (columns) of each nonzero knot span (pages)
    end
    methods
        %% Constructor
        function this = Bspline(varargin)
            % Instantiates a B-spline basis object. Handles polymorphism.
            %  Arguments:
            %   knots: a) number of knot spans; or b) knot span
            %   degree: span-wise polynomial degree of each basis function
            %           (p = 1 will be assummed, if none is given)
            %   smoothness: number of times that each basis function in the
            %               patch can be differentiated (before vanishing, 
            %               i.e. -1 <= smoothness < degree)
            %
            % Default constructor:
            if nargin == 0
                return % leave empty
            end
            % Import knot vector, or generate as specified:
            if length(varargin{1}) == 1 % specified number of non-vanishing knot spans (i.e. CG elements)
                %this.spanCount = varargin{1};
                if varargin{1} > 0
                    this.knots = linspace(0,1,varargin{1}+1);
                else
                    error('Any patch requires, at least, 1 non-vanishing knot span.')
                end
            else
                this.knots = varargin{1}; % given knot vector
                %this.spanCount = length(this.knots) - 1;
            end
            % Set polynomial degree:
            if nargin == 1 || isempty(varargin{2}) || isnan(varargin{2}) || isinf(varargin{2})
                this.degree = 1; % default
            else
                this.degree = varargin{2};
            end
            % Enforce specfied smoothness and/or polynomial degree:
            if this.degree == 0
                this.knots = [-1 1];
            elseif nargin < 3 || isempty(varargin{3}) || isnan(varargin{3})
                this.smoothness = NaN; % undefined (i.e. leave as is)
                this.knots = this.openKnots(this.knots,this.degree);
            else
                [this.knots,this.smoothness] = this.openAndSmoothKnots(...
                    this.knots,this.degree,varargin{3});
            end
            this.knots = this.rescaleKnots(this.knots);
            % Remaining knot vector parameters:
            this.order = this.degree+1;
            this.knotDiffs = diff(this.knots);
            this.knotCount = length(this.knots);
            this.basisCount = this.knotCount-this.order;
            % Connectivity arrays:
            [this.breakCoords,this.span2knot] = unique(this.knots);
            this.span2knot = this.span2knot(2:end)'-1;
            this.nonzeroSpanCount = length(this.span2knot);
            this.span2mode = this.span2knot(1:this.nonzeroSpanCount)' - (this.degree:-1:0);
            % Gauss quadrature data:
            [this.gaussCoords,this.gaussWeights,~] = ...
                Legendre.quadratureGaussLegendre(this.degree);
            % Additional abscissae of interest:
            if this.basisCount == 1
                this.controlCoords = 0;
            else
                this.controlCoords = linspace(-1,1,this.basisCount); % any better idea?
            end
            this.isNodal = (this.degree < 2) || this.isNodal; % p = 0 is Godunov's and p = 1 is nodal FEM with DG coupling
            if this.isNodal
                this.nodeCoords = this.controlCoords';
            end
            % Precomputed basis samples:
            this.left = zeros(this.basisCount,1);
            this.left(1) = 1;
            this.right = flip(this.left);
            % Precomputed operators:
            this.assembleOperators;
            % Sparsity graph:
            this.computeSparsityGraph;
        end
        %% Evaluate basis at locations
        function B = sampleAt(this,x,keepAll)
            % Samples the basis functions at given locations (in parameter
            % space). Vanishing basis functions are removed by default.
            % Arguments
            %  x: 1D row array of evaluation locations
            %  keepAll: (optional, default: false) include vanishing basis
            %           functions in the output
            % Return
            %  l: 2D matrix of evaluated B-spline basis functions (row: basis
            % component; column: evaluation position)
            B = this.samplePatch(this.degree,this.knots,x); % BUG: a sample located at the right-most breakpoint will be sampled at zero
            B = B(:,:,end);
            if nargin == 2 || ~keepAll
%                 toKeep = any(B > 0,2);
%                 toKeep(find(toKeep,1,'first'):find(toKeep,1,'last')) = 1; % BUT WHY?
%                 derivatives = B(toKeep,:);
                B = B(1:end-this.degree,:);
            end
        end
        %% Affine map (from reference patch space to reference knot span space)
        function [sigma,jac] = mapToReferenceKnotSpan(this,xi,spanId)
            % Returns the abcissae in reference knot span space
            % corresponding to the abcissae in reference patch coordinates.
            % 
            % Arguments
            %  xi: abcissae in reference patch coordinates.
            %  spanId: index of the knot span the 'sigma' values are 
            %          referenced to.
            % Output
            %  sigma: abcissae in reference knot span coordinates
            %  jac: Jacobian of the mapping (optional)
            jac = 2/this.knotDiffs(spanId);
            sigma = jac*(xi - this.knots(spanId)) - 1;
        end
        %% Affine map (from reference knot span to reference patch space)
        function [xi,jac] = mapFromReferenceKnotSpan(this,sigma,spanId)
            % Returns the abcissae in reference patch space (-1 <= xi <= 1)
            % corresponding to the abcissae in reference knot span (-1 <=
            % sigma <= 1), for a given knot span (spanId).
            %
            % Arguments
            %  sigma: abcissae in reference knot span coordinates.
            %  spanId: index of the knot span the 'sigma' values are 
            %          referenced to.
            % Output
            %  xi: abcissae in reference patch coordinates
            %  jac: Jacobian of the mapping (optional)
            jac = .5*this.knotDiffs(spanId);
            xi = jac*(sigma+1) + this.knots(spanId);
        end
        %% Compute mass and divergence matrices
        function assembleOperators(this)
            % Computes the components of the mass and gradient matrices
            % of this patch (more or less anologue to local assembly in 
            % FEM). Employs Gauss quadrature at the knot span level, i.e.
            % matrices are exact (to machine accuracy). Additionally,
            % assembles the Vandermonde-like tensors.
            %
            % Outputs
            %  M: Mass matrix, i.e. interior product of each B-spline basis 
            %     function with every other, over the patch domain. 
            %     Row: basis component; column: basis component.
            %  D: Gradient matrix, i.e. interior product over the patch 
            %     domain of each B-spline basis function and every basis 
            %     function's first derivative. Row: basis function;
            %     column: basis function derivative.
            %
            % Preallocation:
            this.vandermonde = zeros(this.order,this.order,this.nonzeroSpanCount);
            this.derivatives = zeros(this.order,this.order,this.nonzeroSpanCount);
            this.massMatrix = sparse(this.basisCount,this.basisCount);
            this.gradientMatrix = this.massMatrix;
            % Loop over non-vanishing spans:
            for l = 1:this.nonzeroSpanCount
                % Left knot of current span:
                n = this.span2knot(l);
                % Affine mapping from reference knot span (sigma) to reference patch (xi):
                [x,spanToElem] = this.mapFromReferenceKnotSpan(this.gaussCoords,n);
                % Vandermonde-like tensors:
                [this.vandermonde(:,:,l),this.derivatives(:,:,l)] = this.sampleSpan(this.knots,this.degree,n,x);
                this.invVandermonde(:,:,l) = inv(this.vandermonde(:,:,l));
                % Knot span contributions to patch-wide inner products:
                j = this.span2mode(l,:); % basis function (global) indices
                for idx = 1:this.order % loop over (local) test function indices
                    r = j(idx); % test function (global) index
                    this.massMatrix(j,r) = (this.vandermonde(:,:,l).*this.vandermonde(idx,:,l))*this.gaussWeights*spanToElem + this.massMatrix(j,r);
                    this.gradientMatrix(j,r) = (this.vandermonde(:,:,l).*this.derivatives(idx,:,l))*this.gaussWeights*spanToElem + this.gradientMatrix(j,r);
                end
            end
        end
        %% Retrieve quadrature points (in reference patch coordinates)
        function x = getGaussCoords(this)
            % NOTE: every knot span has the same Gauss points (in reference
            %       span coordinates), since p is unique throughout the patch.
            x = zeros(1,this.nonzeroSpanCount*(this.degree+1));
            % Loop over knot spans:
            ids = 1:this.degree+1;
            for l = this.span2knot
                x(ids) = this.mapFromReferenceKnotSpan(this.gaussCoords,l);
                ids = ids + this.degree + 1;
            end
        end
    end
    methods (Static)
        %% Re-scale knot vector
        function knots = rescaleKnots(knots,min,max)
            if nargin == 1 % default
                min = -1;
                max = 1;
            end
            knots = (max-min)*(knots-knots(1))/(knots(end)-knots(1))+min;
        end
        %% Isolate interior + edges portion of a given knot vector
        function knots = getInteriorAndEdgeKnots(knots)
            [knots,ids] = unique(knots);
            reps = [diff(ids); 1];
            reps(1) = 1;
            knots = repelem(unique(knots),reps);
        end
        %% Open a given knot vector
        function knots = openKnots(knots,degree)
            knots = DGIGA.getInteriorAndEdgeKnots(knots);
            reps = [degree+1 ones(1,length(knots)-2) degree+1];
            knots = repelem(knots,reps);
        end
        %% Open and smoothen a given knot vector
        function [knots,smoothness] = openAndSmoothKnots(knots,degree,smoothness)
            % Add and/or remove components in a knot vector so that a given
            % smoothness is achieved, ensuring that the result is an open
            % (nonperiodic) knot vector for a given degree.
            knots = unique(knots);
            if isinf(smoothness)
                mult = 1;
            else
                if smoothness > degree-1 && length(knots) > 2
                    smoothness = degree-1;
                    warning(['Enforcing maximum smoothness for a p = ' num2str(degree) ' patch: ' num2str(smoothness) ' times differentiable.'])
                end
                mult = degree-smoothness;
            end
            reps = mult*ones(1,length(knots)-2); % interior knots only
            reps = [degree+1 reps degree+1]; % edge knots also
            knots = repelem(knots,reps);
        end
        %% Consistent L2 projection (norm-preserving)
        function project(mesh,~,fun,q)
            % L2 projection onto a B-spline basis.
            for element = mesh.elements
                basis = element.basis; % handle to current element's basis
                % Employ a projection space of degree q:
                if nargin < 4
                    q = max(50,2*element.basis.degree+1); % default
                end
                [sigmas,w,~] = Legendre.quadratureGaussLegendre(q);
                % Preallocation:
                b = zeros(element.basis.basisCount,length(fun(0)));
                % Assemble vector of initial condition inner products:
                for l = basis.span2knot % loop over non-vanishing knot spans
                    [x,spanToElem] = basis.mapFromReferenceKnotSpan(sigmas,l); % from reference knot span to reference patch coordinates
                    B = basis.sampleSpan(basis.knots,basis.degree,l,x);
                    x = element.mapFromReference(x); % from reference patch to physical patch (i.e. mesh) coordinates
                    f = fun(x');
                    % Inner product between basis functions and the function
                    % being projected, over the range of the current knot
                    % span:
                    r = (l-basis.degree):l; % indices of non-zero basis functions in the current knot span
                    b(r,:) = b(r,:) + B*(w.*f')*spanToElem; % the "elemToMesh" Jacobian cancels out, but "spanToElem" doesn't!
                end
                % Solve for the modal coefficients:
                element.states = b' / element.basis.massMatrix;
            end
        end
        %% Lumped L2 projection (TVD + norm-preserving)
        function projectLumped(mesh,~,fun,q)
            % Lumped L2 projection onto a B-spline basis. From Kuzmin et
            % al, 2010.
            for element = mesh.elements
                basis = element.basis; % handle to current element's basis
                % Employ a projection space of degree q:
                if nargin < 4
                    q = max(50,2*element.basis.degree+1); % default
                end
                [sigmas,w,~] = Legendre.quadratureGaussLegendre(q);
                % Preallocation:
                b = zeros(element.basis.basisCount,length(fun(0)));
                % Assemble vector of initial condition inner products:
                for l = basis.span2knot % loop over non-vanishing knot spans
                    [x,spanToElem] = basis.mapFromReferenceKnotSpan(sigmas,l); % from reference knot span to reference patch coordinates
                    B = basis.sampleSpan(basis.knots,basis.degree,l,x);
                    x = element.mapFromReference(x); % from reference patch to physical patch (i.e. mesh) coordinates
                    f = fun(x');
                    % Inner product between basis functions and the function
                    % being projected, over the range of thre current knot
                    % span:
                    r = (l-basis.degree):l; % indices of non-zero basis functions in the current knot span
                    b(r,:) = b(r,:) + B*(w.*f')*spanToElem; % the "elemToMesh" Jacobian cancels out, but "spanToElem" doesn't!
                end
                % Compute the lumped mass matrix of this element:
                lumpedMass = sum(element.basis.massMatrix,1);
                % Solve for the modal coefficients:
                element.states = b' ./ lumpedMass;
            end
        end
        %% Quasi-interpolatory projection
        function interpolate(mesh,~,fun)
            % Assign control point values from exact function samples. 
            % B-spline function at control locations will, in general, not
            % coincide with the exact one; also, the projection will not be
            % norm-preserving (in general). Thanks to B-spline properties,
            % it will be TVD.
            %
            for element = mesh.elements
                x = element.mapFromReference(element.basis.controlCoords);
                element.states = fun(x);
            end
%             % Apply limiter (if any):
%             if ~isempty(limiter)
%                 limiter.apply(mesh);
%             end
        end
        %% Sample an arbitrary knot vector
        function B = samplePatch(degree,knots,samples)
            % Returns a 3D array of basis function samples (row: knot span; column:
            % sample; page: poly. order). Somewhat vectorized.
            knots = knots(:);
            samples = sort(samples);
            N = length(knots);
            B = zeros(N-1,length(samples),degree+1);
            % Piecewise constants:
            B(:,:,1) = samples >= knots(1:end-1) & samples < knots(2:end);
            ii = find(samples(end) == knots); %%% BUGGY...
            if ~isempty(ii)
                %%%
                id = max(1,ii(1)-1); %%% ¿FIXED?
                B(id,end,1) = 1;
                %%%
                %B(ii(1)-1,end,1) = 1; % special case
            end
            % Higher degrees:
            for p = 1:degree
                j = p+1;
                n = N-j; % index of the right-most non-zero basis function of order j using an (open) N-knot vector
                temp = (samples-knots(1:n))./(knots(j:N-1)-knots(1:n)).*B(1:n,:,p);
                temp(isnan(temp)) = 0;
                B(1:n,:,j) = temp; % ¿ B(1:end-p,:,j) ?
                temp = (knots(1+j:N)-samples)./(knots(1+j:N)-knots(2:n+1)).*B(2:n+1,:,p);
                temp(isnan(temp)) = 0;
                B(1:n,:,j) = B(1:n,:,j) + temp;
            end
        end
        %% Evaluate B-splines and first derivatives on a knot span
        function [vandermonde,derivatives] = sampleSpan(X,p,l,x)
            % Sample the zeroth and first derivatives of all basis 
            % functions of given degree that are non-zero at given 
            % locations (in a given knot span of a given open knot 
            % vector). Based on Piegl & Tiller, 1997: A2.2, and
            % equations 2.5 and 2.7.
            %
            % Arguments
            %  X: knot vector
            %  p: basis functions degree
            %  l: knot span index, such that: X(l) <= x < X(l+1) for all x
            %  x: 1D array of sample locations (column or row)
            % Output
            %  vandermonde: 2D array of all p+1 non-zero basis functions of degree p
            %     (rows), sampled at x (columns)
            %  derivatives: idem, for each function's first derivative
            %
            % Preallocations:
            left = zeros(p,length(x));
            right = left;
            vandermonde = ones(p+1,length(x));
            derivatives = zeros(size(vandermonde));
            % Loop over degrees (p > 0):
            for j=1:p
                left(j,:) = x-X(l+1-j); % left numerator in eq. 2.5
                right(j,:) = X(l+j)-x; % idem, right
                saved0 = zeros(1,length(x));
                saved1 = saved0;
                % Loop over non-zero basis functions of degree j:
                for r=1:j
                    temp = vandermonde(r,:)./(right(r,:)+left(j-r+1,:));
                    vandermonde(r,:) = saved0+right(r,:).*temp;
                    saved0 = left(j-r+1,:).*temp;
                    derivatives(r,:) = saved1-j*temp;
                    saved1 = j*temp;
                end
                vandermonde(j+1,:) = saved0;
                derivatives(j+1,:) = saved1;
            end
        end
        %% Find knot span to which x belongs (open knot vector only)
        function n = findSpan(X,p,x)
            % Returns the knot span to which the given sample point
            % belongs. Not vectorized. Only for open knot vectors.
            % Piegl & Tiller, 1997: A2.1.
            %
            % Arguments:
            %  p: polynomial degree
            %  X: knot vector
            %  x: sample location
            % Output:
            %  n: knot span that contains the sample point x
            %
            % Special case:
            n = length(X)-p-1; % id of the knot at the left of the right breakpoint
            if (x == X(end))
                return;
            end
            % Binary search:
            low = p+1;
            high = n+1;
            n = floor(0.5*(low + high));
            while x < X(n) || x >= X(n+1)
                if x < X(n)
                    high = n;
                else
                    low = n;
                end
                n = floor(0.5*(low + high));
            end
        end
        %% Evaluate non-zero B-splines on a knot span
        function B = basisFuns(X,p,n,x)
            % Sample all basis functions of given degree that are non-zero 
            % at given locations (in a given knot span of a given open knot 
            % vector). Adapted from Piegl & Tiller, 1997: A2.2.
            %
            % Arguments
            %  X: knot vector
            %  p: basis functions degree
            %  n: knot span index, such that: X(n) <= x < X(n+1)
            %  x: array of sample locations
            % Output
            %  B: 2D array of p+1 basis functions of degree p (columns), 
            %     sampled at x (rows); all others (of degree p) are zero
            %
            % Preallocations:
            left = zeros(length(x),p);
            right = left;
            B = ones(length(x),p+1); % trivial case (p = 0)
            % Loop over p > 1 degrees:
            for j=1:p
                left(:,j) = x-X(n+1-j);
                right(:,j) = X(n+j)-x;
                saved = zeros(length(x),1);
                for r=1:j
                    temp = B(:,r)./(right(:,r)+left(:,j-r+1));
                    B(:,r) = saved+right(:,r).*temp;
                    saved = left(:,j-r+1).*temp;
                end
                B(:,j+1) = saved;
            end
        end
        %% Evaluate B-spline basis functions and their derivatives
        function varargout = dersBasisFuns(X,p,n,x)
            % Sample any number of non-zero derivatives (including the 0th 
            % one) at a given location (in a given knot span, of a given
            % open knot vector). Piegl & Tiller, 1997: A2.3.
            %
            % Arguments
            %  X: knot vector
            %  p: basis functions degree
            %  n: knot span index, such that: X(n) <= x < X(n+1)
            %  x: scalar sample location
            % Output
            %  varargout: variable argument list of p+1 column arrays, each 
            %             containing a derivative of p+1 basis functions 
            %             (rows), sampled at the x location.
            %
            % Preallocation:
            ndu = ones(p+1,p+1); % auxiliary array; stores basis function samples (upper triangular part) and knot differences (strictly lower triangular part)
            left = zeros(1,p);
            right = left;
            varargout = cell(1,nargout);
            for j=1:p
                left(j) = x-X(n+1-j);
                right(j) = X(n+j)-x;
                saved = 0;
                for r=1:j
                    % Lower triangular part:
                    ndu(j+1,r) = right(r)+left(j-r+1);
                    temp = ndu(r,j)./ndu(j+1,r);
                    % Upper triangular part:
                    ndu(r,j+1) = saved+right(r).*temp;
                    saved = left(j-r+1).*temp;
                end
                ndu(j+1,j+1) = saved;
            end
            varargout{1} = ndu(:,p+1); % store the 0-th derivatives
            % Compute the remaining derivatives:
            for r=1:p+1 % loop over basis component
                s1 = 1; s2 = 2; % some auxiliary indices
                a = ones(2,p+1); % some auxiliary values
                for k = 1:nargout-1 % loop to compute the k-th derivative
                    d = 0;
                    rk = r-k-1;
                    pk = p-k;
                    if r-1 >= k
                        a(s2,1) = a(s1,1)/ndu(pk+2,rk+1);
                        d = a(s2,1)*ndu(rk+1,pk+1);
                    end
                    if rk >= -1, j1 = 1; else, j1 = -rk; end
                    if r-2 <= pk, j2 = k-1; else, j2 = p-r+1; end
                    for j=j1:j2
                        a(s2,j+1) = (a(s1,j+1) - a(s1,j))/ndu(pk+2,rk+j+1);
                        d = d + a(s2,j+1)*ndu(rk+j+1,pk+1);
                    end
                    if r-1 <= pk
                        a(s2,k+1) = -a(s1,k)/ndu(pk+2,r);
                        d = d + a(s2,k+1)*ndu(r,pk+1);
                    end
                    varargout{k+1}(r,1) = d;
                    j = s1; s1 = s2; s2 = j; % switch rows
                end
            end
            % Multiply by some factors (equation 2.9):
            r = p;
            for k=1:nargout-1
                for j=1:p+1
                    varargout{k+1}(j) = varargout{k+1}(j)*r;
                end
                r = r*(p-k);
            end
        end
    end
end