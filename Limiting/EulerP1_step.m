classdef EulerP1_step < EulerP1
    % Failsafe limiter for Euler equations (high-resolution version).
    %
    % Sets quadratic and above Legendre coefficients (in characteristic
    % variables) to zero, on those elements which experience unphysical
    % states, and "minmods" their slopes.
    %
    % This version acs BOTH after a stage and a step (use with e.g. AFC).
    %
    methods
        %% Constructor
        function this = EulerP1_step(varargin)
            % Superclass constructor:
            this = this@EulerP1(varargin{:});
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