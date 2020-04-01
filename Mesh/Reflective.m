classdef Reflective < Boundary
    % Solid wall boundary that reflects any incoming wave back into the
    % domain. Can be static or moving. See Leveque, 2002, p. 136 and Toro,
    % 2009, p. 496).
    properties (SetAccess = immutable)
        % Speed of the wall, in units of the propagation speed:
        wallSpeed = @(t) 0
    end
    properties (Access = protected)
        % Matrix that, applied to the bound element's degrees 
        % of freedom, gives the corresponding ghost element's DoFs.
        % Encodes a simple reflection (for pressure-like variables).
        A
        % Repetition array that, when applied to a (scalar) wall velocity, 
        % gives a 1D array of velocity DoFs, in the adequate shape to be 
        % added to the ghost states, i.e.: u_ghost = -u_bound*A + u_wall*B
        % (c.f. Toro, 2009, p. 496).
        B
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
            % Computes and stores the LU decomposition of matrix A in
            % A*X = B, where X is the 2D array of ghost states (i.e.
            % Legendre coeffs.), B(1,1:I) is the solution vector at the
            % boundary, and B(2:J,1:I) are the Legendre coefficients of the
            % bound element of this boundary.
            %
            % Set bound and ghost elements (both have the same basis):
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = Element(elements(1).basis,[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                this.ghostElement = Element(elements(end).basis,[elements(end).xR inf]);
            end
            % Set pressure (A) and wall-velocity (B) permutation matrices:
            J = this.boundElement.dofCount;
            if isa(this.ghostElement.basis,'Legendre')
                vals = (-1).^(0:J-1);
                this.A = spdiags(vals',0,J,J);
                this.B = [2 repelem(0,J-1)]; % only the mean value is perturbed
            elseif isa(this.ghostElement.basis,'Lagrange') || isa(this.ghostElement.basis,'BSpline')
                this.A = flip(speye(J));
                this.B = repelem(2,J); % all DoFs are perturbed equally
            else
                error('Basis type not supported.')
            end
        end
        %% Enforce (scalar)
        function apply_scalar(this,physics,solver,~)
            % Computes the ghost element's state vectors matrix by applying
            % the (precomputed) reflection operator. "Polymorphic" on the
            % physics type.
            if isa(physics,'Burgers')
                this.ghostElement.states = -this.ghostElement.states(2,:) + this.wallSpeed(solver.timeNow)*this.B*2; % note: 2*uWall because of 1/2 factor
            elseif isa(physics,'Wave')
                this.ghostElement.states = this.boundElement.states*this.A;
                this.ghostElement.states(2,:) = -this.ghostElement.states(2,:) + this.wallSpeed(solver.timeNow)*this.B;
            elseif isa(physics,'Euler')
                this.ghostElement.states = physics.stateToPrimitive(this.boundElement.states)*this.A; % note: first primitives, then reflection
                this.ghostElement.states(2,:) = -this.ghostElement.states(2,:) + this.wallSpeed(solver.timeNow)*this.B;
                this.ghostElement.states = physics.primitiveToState(this.ghostElement.states);
            else
                error('Physics type not unsupported.')
            end
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (scalar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (u_{wall} = %s)',this.getInfo_scalar@Boundary,func2str(this.wallSpeed));
        end
    end
end