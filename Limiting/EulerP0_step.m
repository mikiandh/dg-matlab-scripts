classdef EulerP0_step < EulerP0
    % Failsafe limiter for Euler equations (Godunov scheme version).
    %
    % Sets linear and above Legendre coefficients (in characteristic
    % variables) to zero, on those elements in which it detects unphysical
    % states.
    %
    % This version acs BOTH after a stage and a step (use with e.g. AFC).
    %
    methods
        %% Constructor
        function this = EulerP0_step(varargin)
            % Superclass constructor:
            this = this@EulerP0(varargin{:});
        end
        %% Apply (override)
        function applyStep(this,mesh,solver)
            this.applyStage(mesh,solver)
        end
        %% Name (extension)
        function name = getName(this)
            name = strrep(this.getName@Limiter,'_step',', also every step');
        end
    end
end