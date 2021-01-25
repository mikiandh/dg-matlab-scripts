classdef EulerP1 < EulerP0
    % Failsafe limiter for Euler equations (high-resolution version).
    %
    % Sets quadratic and above Legendre coefficients (in characteristic
    % variables) to zero, on those elements which experience unphysical
    % states, and "minmods" their slopes.
    %
    methods
        %% Constructor
        function this = EulerP1(varargin)
            % Superclass constructor:
            this = this@EulerP0(varargin{:});
        end
    end
    methods (Access = protected)
        %% Apply equation-wise (override)
        function applyNoSync(this)
            % Limit slopes:
            this.coefs(:,2,:) = permute(this.minmod(...
                permute(this.coefs(:,2,:),[1 3 2]),...
                permute(this.coefs(:,2,:) - this.coefsL(:,2,:),[1 3 2]),...
                permute(this.coefsR(:,2,:) - this.coefs(:,2,:),[1 3 2])...
                ),[1 3 2]);
            % Discard the rest:
            this.coefs(:,3:end,:) = 0;
        end
    end
end