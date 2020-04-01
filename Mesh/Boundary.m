classdef Boundary < handle & matlab.mixin.Heterogeneous
    % Base class from which all valid boundary conditions must derive.
    % Enforcement is made a la weak-Riemann (see Mengaldo et al., 2014), 
    % regardless of the actual type of boundary condition being prescribed.
    properties (SetAccess = protected)
        ghostElement % fictitious (ghost) element
        boundElement % element closest to the domain's bounday, from inside
    end
    methods (Access = protected, Abstract)
        % Initializes the ghost and bound elements of this boundary,
        % according to its type and side (left or right).
        setup_scalar(this,elements,isLeft)
        % Enforces a boundary condition by updating this boundary's ghost
        % element edge values, according to its type and side (left or
        % right).
        apply_scalar(this,physics,solver,isLeft)
    end
    methods (Access = protected)
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
    methods (Static, Access = protected, Sealed)
        function default_object = getDefaultScalarElement
            % Since Boundary is abstract, this functions sets the default
            % object type for heterogeneous Boundary arrays to Periodic.
            default_object = Periodic;
        end
    end
end