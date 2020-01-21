classdef Lagrange < Basis
    properties
        isNodal = true;
        isModal = false;
        isHybrid = false;
        nodeWeights
        vandermonde % nodal_values = modal_coefficients * vandermonde
        invVandermonde % modal_coefficients = nodal_values * invVandermonde
        basisWeights
        derivatives % nodal_derivatives = nodal_coefficients * derivatives
    end
    methods
        %% Constructor
        function lagrange = Lagrange(p)
            if nargin > 0
                lagrange.degree = p;
                lagrange.order = p+1;
                lagrange.basisCount = lagrange.order;
                [lagrange.nodeCoords, lagrange.nodeWeights,lagrange.vandermonde] =...
                    Legendre.quadratureGaussLegendre(p);
                lagrange.gaussCoords = lagrange.nodeCoords;
                lagrange.gaussWeights = lagrange.nodeWeights;
                lagrange.invVandermonde = inv(lagrange.vandermonde);
                lagrange.basisWeights =...
                    lagrange.barycentricWeights(lagrange.nodeCoords);
                lagrange.left = lagrange.sampleAt(-1);
                lagrange.right = lagrange.sampleAt(1);
                lagrange.derivatives = lagrange.derivativeMatrix(...
                    lagrange.nodeCoords,lagrange.basisWeights); % vanilla derivative matrix
                lagrange.assembleMassAndConvectionMatrices;
            end
        end
        %% Evaluate Lagrange polynomial basis
        function l = sampleAt(this,x)
            % See Kopriva, 2009 (algorithm 34, pg. 77).
            % Arguments
            %  x: 1D row array of M evaluation locations
            % Return
            %  l: 2D matrix of evaluated lagrange polynomials (row: basis
            % component; column: evaluation position)
            l = this.basisWeights./((x-this.nodeCoords).*...
                sum(this.basisWeights./(x-this.nodeCoords),1));
            l(x == this.nodeCoords) = 1;
        end
        %% Assemble matrix operators
        function assembleMassAndConvectionMatrices(this)
            % Assembles the mass and convergence matrices in reference
            % element space (i.e. the Jacobian of the mapping between 
            % physical and reference elements is not included in the
            % inner products), for the space spanned by this basis.
            %
            % Preallocation:
            this.massMatrix = diag(this.nodeWeights); % by definition
            this.gradientMatrix = this.derivatives'.*this.nodeWeights; % idem (08/10/2019: it would be practical to redefine 'this.derivatives' as its transposed, but I don't dare change it at this stage)
        end
    end
    methods (Static)
        %% Barycentric weights
        function w = barycentricWeights(x)
            % Adapted from Kopriva, 2009 (algorithm 30, pg. 75).
            % Argument
            %  x: 1D column array of N nodal locations
            % Return
            %  w: 1D column array of N basis component weights
            N = length(x);
            w = 1./prod(eye(N) + repmat(x,1,N) - repmat(x',N,1),2);
        end
        %% Lagrange basis derivative matrix
        function D = derivativeMatrix(x,w)
            % Adapted from Kopriva, 2009 (algorithm 37, pg. 82).
            % Arguments
            %  x: 1D column array of N nodal locations
            %  w: 1D column array of N barycentric weights
            % Return
            %  D: derivative matrix, N X N (row: j-th basis function; 
            %     column: n-th evaluation location)
            %
            % Example: D(1,2) => Derivative of the lagrange basis function
            %          associated to node #1, evaluated at node #2.
            %
            N = length(x);
            % Non-diagonal entries:
            D = w./w'./(x'-x);
            % Diagonal entries (using the 'negative sum trick', pg. 55):
            ids = eye(N,'logical');
            D(ids) = 0;
            D(ids) = -sum(D,1);
        end
        %% Interpolate function on mesh (naive)
        function interpolate(mesh,limiter,fun)
            % Interpolatory projection of the function onto a nodal basis.
            for element = mesh.elements
                x = element.mapFromReference(element.basis.nodeCoords');
                element.states = fun(x);
            end
            % Apply limiter (if any):
            if ~isempty(limiter)
                limiter.apply(mesh);
            end
        end
        %% Project function on mesh (L2)
        function project(mesh,limiter,fun,q)
            % Conservative projection onto a nodal basis.
            for element = mesh.elements
                % Employ a projection space of degree q:
                if nargin < 4
                    q = max(50,2*element.basis.degree+1); % default
                end
                [coords,weights,LGVM] = Legendre.quadratureGaussLegendre(q);
                % Evaluate function at Gauss-Legendre quadrature points:
                x = element.mapFromReference(coords);
                f = fun(x');
                % Project onto a Legendre basis of degree p applying Gauss
                % quadrature with degree q:
                U = (weights'.*f)*...
                    (LGVM(1:element.basis.order,:))';
                U = 0.5*(2*(1:element.basis.order) - 1).*U;
                % Convert from Legendre (modal) to Lagrange (nodal) basis:
                element.states = U*element.basis.vandermonde;
            end
            % Apply limiter (if any):
            if ~isempty(limiter)
                limiter.apply(mesh);
            end
        end
    end
end