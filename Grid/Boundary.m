classdef Boundary < handle & matlab.mixin.Heterogeneous
    % Base class from which all valid boundary conditions must derive.
    % Enforcement is made a la weak-Riemann (see Mengaldo et al., 2014), 
    % regardless of the actual type of boundary condition being prescribed.
    properties (SetAccess = protected)
        ghostElement = Element(DG(0),[-inf inf]) % fictitious (ghost) element
        boundElement % element closest to the domain's bounday, from inside
    end
    methods (Access = protected)
        %% Initialize (scalar, default)
        function setup_scalar(this,elements,isLeft)
            % Assigns bound element and edge coordinate, considering sides.
            if isLeft
                this.boundElement = elements(1);
                this.ghostElement.xR = elements(1).xL;
                this.ghostElement.x = -inf;
            else
                this.boundElement = elements(end);
                this.ghostElement.xL = elements(end).xR;
                this.ghostElement.x = inf;
            end
        end
        %% Enforce (scalar, default)
        function apply_scalar(~,~,~,~)
            % Does nothing.
        end
        %% Information (scalar, default)
        function info = getInfo_scalar(this)
            info = class(this);
        end
    end
    methods (Sealed)
        %% Initialize (vector)
        function setup(these,elements)
            % Saves a handle to the element closest to this boundary and
            % instantiates a ghost element adjacent to it.
            %
            these(1).setup_scalar(elements,true)
            these(2).setup_scalar(elements,false)
        end
        %% Enforce (vector)
        function apply(these,physics,solver)
            % Updates the ghost element of each of these boundaries
            % according to each's scalar apply method.
            these(1).apply_scalar(physics,solver,true)
            these(2).apply_scalar(physics,solver,false)
        end
        %% Information (vector)
        function info = getInfo(these)
            info = sprintf('%s, %s',these(1).getInfo_scalar,these(2).getInfo_scalar);
        end
    end
end