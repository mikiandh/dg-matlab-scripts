classdef Reflective < Boundary
    % Solid wall boundary that reflects any incoming wave back into the
    % domain. Can be static or moving. See Leveque, 2002, p. 136 and Toro,
    % 2009, p. 496).
    properties (SetAccess = immutable)
        % Speed of the wall, in units of the propagation speed:
        wallSpeed = @(t) 0
    end
    methods
        %% Constructor
        function these = Reflective(varargin)
            if nargin > 0
                these(1,nargin) = Reflective;
                [these.wallSpeed] = varargin{:};
                % Convert to (constant) functions of time:
                for this = these
                    if isa(this.wallSpeed,'numeric')
                        this.wallSpeed = str2func(sprintf('@(t) %g',this.wallSpeed));
                    end
                end
            end
        end
    end
    methods (Access = protected)
        %% Initialize (scalar)
        function setup_scalar(this,elements,isLeft)
            % Initializes a p = 0, DG, ghost element.
            %
            % Set bound and ghost elements:
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = Element(DG(0),[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                this.ghostElement = Element(DG(0),[elements(end).xR inf]);
            end
        end
        %% Enforce (scalar)
        function apply_scalar(this,physics,solver,isLeft)
            % Sets the ghost element's edge state vector to match that of
            % the bound element's.
            if isLeft
                state = this.boundElement.stateL;
            else
                state = this.boundElement.stateR;
            end
            if isa(physics,'Burgers')
                state = -state + 2*this.wallSpeed(solver.timeNow);
            elseif isa(physics,'Wave')
                state(2) = -state(2) + 2*this.wallSpeed(solver.timeNow);
            elseif isa(physics,'Euler')
                state = physics.stateToPrimitive(state);
                state(2) = -state(2) + 2*this.wallSpeed(solver.timeNow);
                state = physics.primitiveToState(state);
            else
                error('Physics type not supported.')
            end
            this.ghostElement.states = state;
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (scalar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (u_{wall} = %s)',this.getInfo_scalar@Boundary,func2str(this.wallSpeed));
        end
    end
end