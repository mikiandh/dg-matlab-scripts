classdef ESDGSEM < Lobatto
    methods
        %% Constructor
        function esdgsem = ESDGSEM(varargin)
            esdgsem@Lobatto(varargin{:});
            if nargin > 0
                esdgsem.derivatives = - esdgsem.nodeWeights./...
                    esdgsem.nodeWeights'.*esdgsem.derivatives'; % TO DO
            end
        end
        %% Instantiate from prototype
        function esdgsem = clone(~,p)
            esdgsem = ESDGSEM(p);
        end
    end
    methods (Access = {?Basis,?Element})
        %% DGSEM operator
        function computeResiduals(this,element,physics) % TO DO
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