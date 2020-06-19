classdef DGIGA_nodal_AFC < DGIGA_AFC
    methods
        %% Constructor
        function this = DGIGA_nodal_AFC(varargin)
            this@DGIGA_AFC(varargin{:});
        end
        %% Instantiate from prototype
        function this = clone(prototype,degree)
            this = DGIGA_nodal_AFC(prototype.knots,degree,prototype.smoothness);
        end
    end
    methods (Access = {?Basis,?Element})
        %% DGIGA-AFC operator (low-order predictor)
        function computeResiduals(this,element,physics)
            element.residuals = element.states;
            element.fluxes = physics.flux(element.residuals);
            element.residuals = ...
                element.fluxes*this.gradientMatrix...
                - element.riemannR.*this.right'...
                - element.riemannL.*this.left';
            this.diffuseResiduals(element,physics)
            element.residuals = 2/element.dx*element.residuals./this.lumpedMassMatrixDiagonal;
        end
    end
end