classdef SSP_RK3 < TimeIntegrator
    properties (Constant)
        order = 3;
        stageCount = 3;
        amplificationFactorFun = @(z) 1 + z + 1/2*z.^2 + 1/6*z.^3;
    end
    methods
        %% Constructor
        function ssp_rk3 = SSP_RK3(varargin)
            ssp_rk3@TimeIntegrator(varargin{:});
        end
        %% RK stage update
        function applyStage(this,element,stage)
            switch stage
                case 1
                    element.extraStates = element.states;
                    element.states = element.states + this.timeDelta * element.residuals;
                case 2
                    element.states = 0.75*element.extraStates + 0.25*(element.states + this.timeDelta*element.residuals);
                otherwise
                    element.states = 1/3*element.extraStates + 2/3*(element.states + this.timeDelta*element.residuals);
            end
        end
    end
end