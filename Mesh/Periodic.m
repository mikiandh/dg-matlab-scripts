classdef Periodic < Boundary
    % Periodic coupling of this boundary with that at the opposite end of 
    % the mesh.
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
                this.ghostElement = elements(end);
            else
                this.boundElement = elements(end);
                this.ghostElement = elements(1);
            end
        end
        %% Enforce (scalar)
        function apply_scalar(varargin)
            % Intentionally does nothing (the ghost element is a handle to
            % an actual element, which is updated by the solver).
        end
    end
end