classdef DG < Legendre
    properties
        norms % column array containing the (inverse of) non-zero interior products between basis components per unit length
    end
    methods
        %% Constructor
        function dg = DG(varargin)
            dg@Legendre(varargin{:});
            if nargin > 0
                dg.norms = 2*(1:dg.order)-1;
                dg.derivatives =...
                    dg.gaussWeights.*dg.derivatives'; % modified derivative matrix (weighted 1st derivative, row: basis comp.; column: Gauss point)
            end
        end
        %% Instantiate from prototype
        function dg = clone(~,p)
            dg = DG(p);
        end
        %% DG operator
        function computeResiduals(this,element,physics)
            element.residuals = element.states*this.vandermonde;
            element.residuals = physics.flux(element.residuals);
            element.residuals = - element.residuals*this.derivatives;      % volume contribution
            element.residuals = element.residuals +...                     % edge contributions
                element.riemannR.*this.right' +...
                element.riemannL.*this.left';
            element.residuals = - this.norms./element.dx.*element.residuals;
        end
    end
    methods (Static)
        %% Interpolate function onto mesh
        function interpolate(mesh,limiter,fun)
            % Interpolatory projection of the function onto a MODAL basis.
            for element = mesh.elements
                x = element.mapFromReference(element.basis.gaussCoords');
                element.states = fun(x);
                element.nodalToModal;
            end
            % Apply limiter (if any):
            if ~isempty(limiter)
                limiter.apply(mesh);
            end
        end
    end
    
end