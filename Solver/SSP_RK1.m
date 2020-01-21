classdef SSP_RK1 < TimeIntegrator
    properties (Constant)
        order = 1;
        stageCount = 1;
        amplificationFactorFun = @(z) 1+z;
    end
    methods
        %% Constructor
        function ssp_rk1 = SSP_RK1(varargin)
            ssp_rk1@TimeIntegrator(varargin{:});
        end
        %% RK stage update
        function applyStage(this,element,~)
            element.states = element.states + this.timeDelta * element.residuals;
        end
    end
end