classdef DGIGA < Bspline
    methods
        %% Constructor
        function this = DGIGA(varargin)
            this@Bspline(varargin{:});
        end
        %% Instantiate from prototype
        function this = clone(prototype,degree)
            this = DGIGA(prototype.knots,degree,prototype.smoothness);
        end
    end
    methods (Access = {?Basis,?Element})
        %% DGIGA operator
        function computeResiduals(this,element,physics)
            element.residuals = element.states*this.controlVandermonde;
            element.fluxes = physics.flux(element.residuals)/this.controlVandermonde;
            element.residuals = (element.fluxes*this.gradientMatrix...
                - element.riemannR.*this.right'...
                - element.riemannL.*this.left');
            element.residuals = element.residuals / this.massMatrix;
            element.residuals = element.residuals*2/element.dx; % jacobian of the mapping
        end
    end
end