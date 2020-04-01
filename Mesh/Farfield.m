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
        %% Initialize (scalar)
        function setup_scalar(this,elements,isLeft)
            % Sets this boundary's bound and ghost elements. The latter
            % uses a FV basis (i.e. element-wide constant solution, p = 0).
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = Element(DG(0),[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                this.ghostElement = Element(DG(0),[elements(end).xR inf]);
            end
        end
        %% Enforce (scalar)
        function apply_scalar(this,~,solver,~)
            % Updates the ghost state vector according to this boundary's 
            % farfield state function (of time). Propagates it to this 
            % boundary's ghost element's edges.
            this.ghostElement.states = this.states(solver.timeNow);
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (sclar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (q_{\\infty} = %s)',this.getInfo_scalar@Boundary,func2str(this.states));
        end
    end
end