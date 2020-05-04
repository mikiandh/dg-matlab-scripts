classdef Transmissive < Boundary
    % Fictitious boundary, i.e. allows any characteristic to leave without
    % reflecting any spurious wave back into the domain (see Leveque, 2002,
    % p. 134 and Toro, 2009, p. 495). Also known as: zero-gradient, 
    % non-reflective, absorbing.
    properties (Access = protected)
        % Matrix that couples the ghost and bound element's state vectors 
        % to each other: Q_ghost = Q_bound*C
        C
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
        %% Initialize (scalar)
        function setup_scalar(this,elements,isLeft)
            % Computes the constraint matrix.
            %
            % Set bound and ghost elements (both have the same #DoFs):
            if isLeft
                this.boundElement = elements(1);
                N = this.boundElement.dofCount;
                p = N-1;
                s = [1 -1];
                this.ghostElement = Element(DG(p),[-inf elements(1).xL]);
            else
                this.boundElement = elements(end);
                N = this.boundElement.dofCount;
                p = N-1;
                s = [-1 1];
                this.ghostElement = Element(DG(p),[elements(end).xR inf]);
            end
            % Assemble the matrix:
            AB = permute(Legendre.getDerivatives(N,s,p),[1 3 2]);
            this.C = AB(:,:,2)/AB(:,:,1);
        end
        %% Enforce (scalar)
        function apply_scalar(this,~,~,~)
            % Updates the ghost element's state vectors such that the
            % approximate solution is smooth at the boundary.
            %
            % Right-hand side matrix:
            this.ghostElement.states = this.boundElement.getLegendre*this.C;
            % Evaluate ghost element at its edges:
            this.ghostElement.interpolateStateAtEdges
        end
    end
end