classdef Krivodonova < BDF
    % Hierarchical moment limiter by Krivodonova (2007). Relaxes the
    % minmod limiter by using less restrictive approximated higher-order 
    % moments.
    %
    methods
        %% Constructor
        function this = Krivodonova(varargin)
            % Superclass constructor:
            this = this@BDF(varargin{:});
        end
    end
    methods (Access = protected)
        %% Apply equation-wise (override)
        function applyNoSync(this)
            % Applies the minmod limiter asynchronously (i.e. on each
            % characteristic variable separately).
            %
            % Loop over characteristic variables:
            for i = 1:this.I
                % Flag every element for limiting:
                k = true(this.K,1);
                % Limit hierarchically (top to bottom):
                for j = this.J:-1:1
                    % Get limited coefs.:
                    coefs_min = this.coefs(k,j+1,i);
                    this.coefs(k,j+1,i) = this.minmod(this.coefs(k,j+1,i),...
                        this.coefs(k,j,i) - this.coefsL(k,j,i),...
                        this.coefsR(k,j,i) - this.coefs(k,j,i));
                    % Update list of elements flagged for limiting:
                    k(k) = abs(this.coefs(k,j+1,i) - coefs_min) > 1e-10 | ~coefs_min;
                end
            end
        end
    end
end