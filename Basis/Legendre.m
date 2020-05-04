classdef Legendre < Basis
    properties (Constant)
        isNodal = false
        isModal = true
        isHybrid = false
    end
    properties
        vandermonde % matrix of basis components (rows) sampled at p+1 Gauss quadrature points (columns)
        invVandermonde % matrix of Gauss quadrature point contributions (rows) to modal components (columns)
        derivatives % matrix of first derivative of basis components (rows) sampled at p+1 Gauss quadrature points (columns)
    end
    methods
        %% Constructor
        function this = Legendre(p)
            if nargin > 0
                this.degree = p;
                this.order = p+1;
                this.basisCount = this.order;
                this.breakCoords = [-1 1];
                [this.gaussCoords, this.gaussWeights,this.vandermonde] = this.quadratureGaussLegendre(p);
                this.dofCoords = this.gaussCoords;
                this.invVandermonde = inv(this.vandermonde);
                [this.left,this.right] = ...
                    this.getComponentsAtEdges(this.degree);
                this.derivatives = ...
                    this.diffLegendre(this.gaussCoords',this.vandermonde);
                this.assembleMassAndConvectionMatrices;
                % Sparsity graph:
                this.computeSparsityGraph;
            end
        end
        %% Evaluate Legendre polynomial basis
        function l = sampleAt(this,x)
            % Samples the basis at given locations (in reference space).
            % Argument
            %  x: 1D row array of evaluation locations
            % Return
            %  l: 2D matrix of evaluated legendre polynomials (row: basis
            % component; column: evaluation position)
            [~,l] = this.getLegendreAndDerivatives(this.degree+1,x);
        end
        %% Assemble matrix operators
        function assembleMassAndConvectionMatrices(this)
            % Assembles the mass and convergence matrices in reference
            % element space (i.e. the Jacobian of the mapping between 
            % physical and reference elements is not included in the
            % inner products), for the space spanned by this basis.
            %
            % Preallocation:
            this.massMatrix = zeros(this.basisCount);
            this.gradientMatrix = this.massMatrix;
            for r = 1:this.basisCount
                this.massMatrix(:,r) = (this.vandermonde.*this.vandermonde(r,:))*this.gaussWeights;
                this.gradientMatrix(:,r) = (this.vandermonde.*this.derivatives(r,:))*this.gaussWeights;
            end
        end
        %% Interpolatory projection
        function interpolate(this,element,fun0)
            % Interpolatory projection of the function onto the nodal
            % (i.e. Lagrange) counterpart of this basis.
            %
            element.states = fun0(element.getGaussCoords')*this.invVandermonde;
        end
    end
    methods (Static)
        %% Legendre polynomials and 1st derivative at given points
        function [dPn,Pn] = getLegendreAndDerivatives(n,xi)
            % Returns evaluations of the Legendre polynomials and their
            % first derivative, for degrees from 0 to n-1 at xi. See
            % Kopriva, algorithm 22 (pp. 63).
            %
            % Arguments
            %  n: highest-degree polynomial to sample (n = 1 <-> zero degree)
            %  xi: evaluation locations (1D column array)
            %
            % Return
            %  dPn: derivative of Legendre polynomials 1 to n (n x length(xi))
            %  Pn: Legendre polynomial 1 to n
            %
            % Preallocate:
            dPn = zeros(n,length(xi));
            Pn = ones(n,length(xi));
            % Trivial cases:
            if n == 1
                return;
            end
            dPn(2,:) = 1;
            Pn(2,:) = xi;
            % Sample the Legendre polynomials and their derivatives:
            for j = 2:n-1
                Pn(j+1,:) = ((2*j-1)*xi.*Pn(j,:) - (j-1)*Pn(j-1,:))/j;
                dPn(j+1,:) = j*Pn(j,:) + xi.*dPn(j,:);
            end
        end
        %% Legendre polynomial derivatives (known Vandermonde matrix)
        function dPn = diffLegendre(xi,Pn)
            % P'_0:
            dPn = zeros(size(Pn)); if size(Pn,1) == 1; return; end
            % P'_k, k > 0:
            dPn(2,:) = 1;
            for j = 2:size(Pn,1)-1
                dPn(j+1,:) = j*Pn(j,:) + xi.*dPn(j,:);
            end
        end
        %% Derivatives at given points
        function p = getDerivatives(n,xi,s)
            % Returns evaluations of 1 to s derivatives of Legendre
            % polynomials 1 to n, at a set of points.
            % Based on Kopriva, 2009 (p. 24).
            %
            % Arguments
            %  n: highest-degree polynomial to sample (n = 1 <-> zero degree)
            %  xi: evaluation locations (row array)
            %  s: highest order derivative to compute (s = 0 <-> polynomial itself)
            %
            % Return
            %  p: 0 to s derivatives (pages) of 1 to n Legendre polynomials (rows) at points xi (columns)
            %
            % Preallocate:
            p = zeros(n,length(xi),s+1);
            % Trivial, nonzero cases:
            p(1,:,1) = 1;
            if n > 1
                p(2,:,1) = xi;
                if s > 0
                    p(2,:,2) = 1;
                end
            end
            % Zeroth derivative (i.e. polynomials themselves):
            for j = 2:n-1
                p(j+1,:,1) = ((2*j-1)*xi.*p(j,:,1) - (j-1)*p(j-1,:,1))/j; % (j+1)th degree Legendre polynomial
            end
            % Higher-order derivatives:
            for k = 1:s
                for j = 2:n-1
                    % Differentiated Bonnet's formula:
                    p(j+1,:,k+1) = (2*j-1)*p(j,:,k) + p(j-1,:,k+1); % kth derivative of the (j+1)th degree Legendre polynomial
                end
            end
        end
        %% Legendre basis sampled at edges
        function [left, right] = getComponentsAtEdges(degree)
            % Argument
            %  degree: degree of the Legendre basis to sample
            % Outputs
            %  left, right: column array (row: basis component) of samples 
            %               at each edge (left or right)
            left = (-1).^(0:degree)';
            right = ones(degree+1,1);
        end
        %% Legendre-Gauss quadrature data
        function [coords,weights,LGVM] = quadratureGaussLegendre(N)
            % Arguments
            %  N: Gauss quadrature degree
            % Returns
            %  coords: position of N+1 Gauss-Legendre nodes (left to right)
            %  weights: N+1 Gauss-Legendre weights
            %  LGVM: Gauss-Legendre vandermonde matrix, i.e. modal to
            %        nodal operator (N+1 legendre comp. x N+1 gauss nodes)
            if N == 0
                coords = 0;
                weights = 2;
                LGVM = 1;
                return;
            end
            N1=N+1;
            N2=N+2;
            xu=linspace(-1,1,N1)';
            % Initial guesses:
            coords=-cos((2*(0:N)'+1)*pi/(2*N+2))-(0.27/N1)*sin(pi*xu*N/N2);
            LGVM=zeros(N1,N2); % Legendre-Gauss Vandermonde Matrix
            dLGV=zeros(N1,N2); % derivative of LGVM
            % Compute the zeros of the N-th Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method:
            y0=2;
            while max(abs(coords-y0))>eps
                LGVM(:,1)=1;
                LGVM(:,2)=coords;
                for k=2:N1
                    LGVM(:,k+1)=( (2*k-1)*coords.*LGVM(:,k)-(k-1)*LGVM(:,k-1) )/k;
                end
                dLGV=(N2)*( LGVM(:,N1)-coords.*LGVM(:,N2) )./(1-coords.^2);
                y0=coords;
                coords=y0-LGVM(:,N2)./dLGV;
            end
            weights = 2./((1-coords.^2).*dLGV.^2)*(N2/N1)^2; % compute weights
            LGVM = (LGVM(:,1:end-1))'; % clean-up the Vandermonde matrix
        end
    end
    methods (Access = {?Basis,?Element})
        %% Identity projection (forwards)
        function modes = getLegendre(~,element,j,i)
            % Returns selected modes of the given element.
            switch nargin               
                case 2
                    modes = element.states;
                case 3
                    modes = element.states(:,j);
                case 4
                    modes = element.states(i,j);
            end
        end
        %% Identity projection (backwards)
        function setLegendre(~,element,modes,i)
            % Sets given Legendre modes into the given element.
            switch nargin
                case 3
                    element.states = modes;
                case 4
                    element.states(i,:) = modes;
            end
        end
        %% Legendre to Lagrange projection
        function nodes = getLagrange(this,element,j,i)
            % Returns the expansion coefficients of the projection of this
            % basis to a Lagrange one of the same length.
            switch nargin
                case 2
                    nodes = element.states*this.vandermonde;
                case 3
                    nodes = element.states*this.vandermonde(:,j);
                case 4
                    nodes = element.states(i,:)*this.vandermonde(:,j);
            end
        end
        %% Lagrange to Legendre projection
        function setLagrange(this,element,nodes,i)
            % Sets given Lagrange expansion coefficients to selected
            % state array entries; assumes that the Lagrange and the 
            % element's bases have the same length.
            switch nargin
                case 3
                    element.states = nodes*this.invVandermonde;
                case 4
                    element.states(i,:) = nodes*this.invVandermonde;
            end
        end
    end
end