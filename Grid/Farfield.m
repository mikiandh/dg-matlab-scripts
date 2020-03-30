classdef Farfield < Boundary
    % Specified inflow or outflow (no distinction) state. Weak-Riemann
    % prescription guarantees well-posedness (see Mengualdo et al., 2014,
    % page 7).
    properties (SetAccess = immutable)
        states = @(t) 0
    end
    methods
        %% Constructor
        function these = Farfield(varargin)
            if nargin > 0
                these(1,nargin) = Farfield;
                [these.states] = varargin{:};
                % Convert to (constant) functions of time:
                for this = these
                    if isa(this.states,'numeric')
                        this.states = str2func(['@(t) [' num2str(this.states(:)') ']'' ']);
                    end
                end
            end
        end
    end
    methods (Access = protected)
        %% Enforce (scalar, extension)
        function apply_scalar(this,~,solver,~)
            this.ghostElement.states = this.states(solver.timeNow);
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (sclar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (q_{Inf} = %s)',this.getInfo_scalar@Boundary,func2str(this.states));
        end
    end
end