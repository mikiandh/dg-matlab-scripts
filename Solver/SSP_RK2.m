classdef SSP_RK2 < Solver
    properties (Constant)
        order = 2;
        stageCount = 2;
        amplificationFactorFun = @(z) 1 + z + 1/2*z.^2;
    end
    methods
        %% Constructor
        function ssp_rk2 = SSP_RK2(varargin)
            ssp_rk2@Solver(varargin{:});
        end
        %% RK stage update
        function applyStage(this,element,stage)
            switch stage
                case 1
                    element.extraStates = element.states;
                    element.states = element.states + this.timeDelta * element.residuals;
                otherwise
                    element.states = .5*(element.extraStates + element.states + this.timeDelta*element.residuals);
            end
        end
    end
end