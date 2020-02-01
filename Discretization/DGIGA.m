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
            element.residuals = element.residuals / this.massMatrix;
%%% TEST 1 (bad arrangement)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             TEST1 = vertcat(element.states(1,:)',element.states(2,:)');
%             K = [zeros(this.basisCount) 1^2*this.gradientMatrix'
%                  this.gradientMatrix'   zeros(this.basisCount)];
%             TEST1 = K*TEST1;
%             TEST1 = (reshape(TEST1',this.basisCount,2))';
%             TEST1 == element.fluxes*this.gradientMatrix
%%% TEST 2 (good arrangement)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             TEST2 = element.states(:);
%             K = zeros(2*this.basisCount);
%             for r = 1:this.basisCount
%                 for j = 1:this.basisCount
%                     rs = 2*(r-1)+(1:2);
%                     js = 2*(j-1)+(1:2);
%                     K(rs,js) = [0 1^2; 1 0]*this.gradientMatrix(j,r);
%                 end
%             end
%             TEST2 = K*TEST2;
%             TEST2 = reshape(TEST2,2,this.basisCount);
%             TEST2 == element.fluxes*this.gradientMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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