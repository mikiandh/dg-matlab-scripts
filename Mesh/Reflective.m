classdef Reflective < Boundary
    % Solid wall boundary that reflects any incoming wave back into the
    % domain. Can be static or moving. See Leveque, 2002, p. 136 and Toro,
    % 2009, p. 496). Smoothness-preserving treatment.
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
            % Sets bound and ghost elements (both have the same basis), as
            % well as this boundary's auxiliary matrices A, B and C.
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = Element(elements(1).basis,[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                this.ghostElement = Element(elements(end).basis,[elements(end).xR inf]);
            end
        end
        %% Enforce (scalar)
        function apply_scalar(this,physics,solver,~)
            % Computes the ghost element's state vectors matrix by applying
            % a discrete reflection operator. "Polymorphic" on the physics
            % type.
            %
            % Physics-dependent constraint:
            switch class(physics)
                case 'Burgers'
                    this.ghostElement.states = -this.boundElement.getLagrange + 4*this.wallSpeed(solver.timeNow);
                case 'Wave'
                    this.ghostElement.states = this.boundElement.getLagrange;
                    this.ghostElement.states(2,:) = -this.ghostElement.states(2,:) + 2*this.wallSpeed(solver.timeNow);
                case 'Euler'
                    this.ghostElement.states = Euler.stateToPrimitive(this.boundElement.getLagrange);
                    this.ghostElement.states(2,:) = -this.ghostElement.states(2,:) + 2*this.wallSpeed(solver.timeNow);
                    this.ghostElement.states = Euler.primitiveToState(this.ghostElement.states);
                otherwise
                    error('Physics not supported.')
            end
            % Reflect (as nodal) and set into ghost element (as modal):
            this.ghostElement.setLagrange(this.ghostElement.states(:,end:-1:1));
            % Evaluate at edges:
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (scalar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (u_{wall} = %s)',this.getInfo_scalar@Boundary,func2str(this.wallSpeed));
        end
    end
end