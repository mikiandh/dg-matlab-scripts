classdef DGSEM < Lagrange
    methods
        %% Constructor
        function dgsem = DGSEM(varargin)
            dgsem@Lagrange(varargin{:});
            if nargin > 0
                dgsem.derivatives = - dgsem.nodeWeights./...
                    dgsem.nodeWeights'.*dgsem.derivatives'; % modified derivative matrix
                dgsem.assembleFourierMatrices;
            end
        end
        %% Instantiate from prototype
        function dgsem = clone(~,p)
            dgsem = DGSEM(p);
        end
    end
    methods (Access = {?Basis,?Element})
        %% DGSEM operator
        function computeResiduals(this,element,physics)
            element.computeFluxesFromStates(physics);
            element.residuals = element.fluxes*this.derivatives; % volume contribution
            element.residuals = element.residuals +...                         % edge contributions
                (element.riemannR.*this.right' +...
                element.riemannL.*this.left')./...
                this.nodeWeights';
            element.residuals = - 2/element.dx*element.residuals;
        end
    end
end