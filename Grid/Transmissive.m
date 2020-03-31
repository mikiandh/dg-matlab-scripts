classdef Transmissive < Boundary
    % Fictitious boundary, i.e. allows any characteristic to leave without
    % reflecting any spurious wave back into the domain (see Leveque, 2002,
    % p. 134 and Toro, 2009, p. 495). Also known as: zero-gradient, 
    % non-reflective, absorbing.
    properties (Access = protected)
        % Decomposition A = P'*L*U, matrix which couples the ghost
        % element's state vectors with the bound element Legendre 
        % coefficients and left/right edge value, so that the zero-gradient
        % constraint is satisfied (more or less exactly) for all 
        % derivatives of the solution at the boundary.
        L,U,P
    end
    methods
        %% Constructor
        function these = Transmissive(n)
            if nargin > 0
                these(1,n) = Transmissive;
            end
        end
    end
    methods(Access = protected)
        %% Initialize (scalar, extension)
        function setup_scalar(this,elements,isLeft)
            % Computes and stores the LU decomposition of matrix A in
            % A*X = B, where X is the 2D array of ghost states (i.e.
            % Legendre coeffs.), B(1,1:I) is the solution vector at the
            % boundary, and B(2:J,1:I) are the Legendre coefficients of the
            % bound element of this boundary.
            %
            % Base class setup:
            this.setup_scalar@Boundary(elements,isLeft)
            % Sparse constraint matrix assembly:
            J = this.boundElement.dofCount; % number of rows/columns
            s = isLeft - ~isLeft; % +1 if left boundary, -1 if right
            A = sparse([repelem(1,J) 2:J],[1:J 2:J],[s.^(2:J+1) repelem(1,J-1)]);
            % LU decomposition:
            [this.L,this.U,this.P] = lu(A);
        end
        %% Enforce (scalar, override)
        function apply_scalar(this,~,~,isLeft)
            % Updates the ghost element's state vectors such that the
            % approximate solution is smooth at the boundary. Employs an
            % LU-factorized constraint matrix, such that: x = U\(L\P*b).
            %
            % Right-hand side matrix:
            this.ghostElement.states = this.boundElement.getLegendre;
            if isLeft
                this.ghostElement.states(:,1) = this.boundElement.stateL;
            else
                this.ghostElement.states(:,1) = this.boundElement.stateR;
            end
            % Solve:
            this.ghostElement.states = (this.U\(this.L\this.P*this.ghostElement.states'))';
            % Evaluate ghost element at its edges:
            this.ghostElement.interpolateStateAtEdges
        end
    end
end