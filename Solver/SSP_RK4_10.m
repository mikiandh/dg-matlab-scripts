classdef SSP_RK4_10 < Solver
    properties (Constant)
        order = 4;
        stageCount = 10;
        amplificationFactorFun = @(z) z.^10/251942400 + z.^9/4199040 + z.^8/155520 + z.^7/9720 + (7*z.^6)/6480 + (17*z.^5)/2160 + z.^4/24 + z.^3/6 + z.^2/2 + z + 1;
    end
    methods
        %% Constructor
        function ssp_rk4_10 = SSP_RK4_10(varargin)
            ssp_rk4_10@Solver(varargin{:});
        end
        %% RK stage update
        function applyStage(this,element)
            switch this.stageNow
                case 1
                    element.extraStates{2} = element.states;
                    element.states = element.states + 1/6*this.timeDelta * element.residuals;
                case 5
                    element.extraStates{1} = element.states; % state after the 4th stage
                    element.extraResiduals = element.residuals; % residual of the state after the 4th stage
                    element.states = 0.6*element.extraStates{2} + 0.4*element.states + 1/15*this.timeDelta*element.residuals;
                case 10
                    element.states = 0.04*element.extraStates{2} + 0.36*element.extraStates{1} + 0.6*element.states + 0.06*this.timeDelta*element.extraResiduals + 0.1*this.timeDelta*element.residuals;
                otherwise
                    element.states = element.states + 1/6*this.timeDelta * element.residuals;
            end
        end
    end
end