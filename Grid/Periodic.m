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
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement = elements(end);
            else
                this.boundElement = elements(end);
                this.ghostElement = elements(1);
            end
        end
    end
end