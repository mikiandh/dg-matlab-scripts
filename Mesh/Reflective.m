classdef Reflective < Boundary
    % Solid wall boundary that reflects any incoming wave back into the
    % domain. Can be static or moving. See Leveque, 2002, p. 136 and Toro,
    % 2009, p. 496).
    properties (SetAccess = immutable)
        wallSpeed = @(t) 0 % nondim. by the propagation speed (i.e. Mach number)
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
            % Set bound and ghost elements (both have the same #DoFs):
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = Element(DG(elements(1).dofCount-1),[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                this.ghostElement = Element(DG(elements(end).dofCount-1),[elements(end).xR inf]);
            end
            % Sparse constraint matrix assembly:
            J = this.boundElement.dofCount; % number of rows/columns
            s = isLeft - ~isLeft; % +1 if left boundary, -1 if right
            A = sparse([repelem(1,J) 2:J],[1:J 2:J],[s.^(2:J+1) repelem(1,J-1)]);
            % LU factorization:
            [this.L,this.U,this.P] = lu(A);
        end
        %% Enforce (scalar)
        function apply_scalar(this,physics,solver,isLeft)
            if isLeft
                state = this.boundElement.stateL;
            else
                state = this.boundElement.stateR;
            end
            this.ghostElement.states = Reflective.(class(physics))(state,this.wallSpeed(solver.timeNow));
            this.ghostElement.interpolateStateAtEdges
        end
        %% Information (scalar, extension)
        function info = getInfo_scalar(this)
            info = sprintf('%s (u_{wall} = %s)',this.getInfo_scalar@Boundary,func2str(this.wallSpeed));
        end
    end
    methods (Static, Access = protected)
        %% Ghost state (advection)
        function q = Advection(~,~) %#ok<STOUT>
            error('Reflective boundary condition is not defined for the Advection equation.')
        end
        %% Ghost state (Burgers)
        function q = Burgers(q,uWall)
            q = -q + 4*uWall;
        end
        %% Ghost state (Wave)
        function q = Wave(q,uWall)
            q(2) = -q(2) + 2*uWall;
        end
        %% Ghost state (Euler)
        function q = Euler(q,uWall)
            q = Euler.stateToPrimitive(q);
            q(2) = -q(2) + 2*uWall;
            q = Euler.primitiveToState(q);
        end
    end
end