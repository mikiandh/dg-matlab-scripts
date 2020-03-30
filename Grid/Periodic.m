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
        %% Initialize (scalar, extension)
        function setup_scalar(this,elements,isLeft)
            % Replaces the ghost element adjacent to a boundary with the
            % element adjacent to the opposite one.
            %
            this.setup_scalar@Boundary(elements,isLeft);
            if isLeft
                this.ghostElement = elements(end);
            else
                this.ghostElement = elements(1);
            end
        end
    end
end