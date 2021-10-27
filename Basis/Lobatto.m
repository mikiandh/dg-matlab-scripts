classdef Lobatto < Lagrange
    methods
        %% Constructor
        function this = Lobatto(varargin)
            if nargin > 0
                if varargin{1} < 1
                    error('Gauss-Lobatto nodes are not defined for p < 1.')
                end
            end
            this@Lagrange(varargin{:});
        end
        %% Initialize Gauss-Lobatto quadrature variables
        function setQuadratureVars(this)
            [this.nodeCoords, this.nodeWeights, this.vandermonde] =...
                Legendre.quadratureGaussLobatto(this.degree);
        end
    end
end