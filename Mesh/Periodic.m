classdef Periodic < Boundary
    % Periodic coupling of this boundary with that at the opposite end of 
    % the mesh. Done by storing both a (deep) copy of the bound element of
    % the opposite boundary, and a handle to it (for updating purpouses).
    properties (Access = protected)
        % Handle to the "mesh-owned" instance of this boundary's ghost
        % element:
        oppositeBoundElement
    end
    methods
        %% Constructor
        function these = Periodic(n)
            if nargin > 0
                these(1,n) = Periodic;
            end
        end
    end
    methods (Access = protected)
        %% Initialize (scalar)
        function setup_scalar(this,elements,isLeft)
            % Sets the ghost element adjacent to this boundary to a handle
            % of the element adjacent to the domain's opposite boundary.
            if isLeft
                this.boundElement = elements(1);
                this.oppositeBoundElement = elements(end);
                this.ghostElement = Element(elements(end).basis,[elements(end).xL elements(end).xR]);
            else
                this.boundElement = elements(end);
                this.oppositeBoundElement = elements(1);
                this.ghostElement = Element(elements(1).basis,[elements(1).xL elements(1).xR]);
            end
        end
        %% Enforce (scalar)
        function apply_scalar(this,varargin)
            % Updates the ghost element's states with those of its "real"
            % counterpart (i.e. the copy of it owned by the mesh).
            this.ghostElement.states = this.oppositeBoundElement.states;
            this.ghostElement.interpolateStateAtEdges
        end
    end
end