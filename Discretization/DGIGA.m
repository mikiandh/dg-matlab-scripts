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
        %% DGIGA operator
        function computeResiduals(this,element,physics)
            element.computeFluxesFromStates(physics);
            element.residuals = 2/element.dx*(...
                element.fluxes*this.gradientMatrix - ...
                element.riemannR.*this.right' - ...
                element.riemannL.*this.left');
            %element.residuals = element.residuals*this.invMassMatrix;
            element.residuals = element.residuals / this.massMatrix;
        end
    end
    methods (Static)
        %% Constrained projection
        function projectAFC(varargin)
        % Uses the constrained L2 projection from DGIGA_AFC.
        DGIGA_AFC.project(varargin{:});
        end
        %% Constrained projection (P, Q coupling)
        function projectAFC_Matthias(varargin)
        % Uses the constrained L2 projection from DGIGA_AFC that employs
        % Matthias's coupling proposal.
        DGIGA_AFC.project_Matthias(varargin{:});
        end
    end
end