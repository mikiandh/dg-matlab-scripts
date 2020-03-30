classdef Transmissive < Boundary
    % Fictitious boundary, i.e. allows any characteristic to leave without
    % reflecting any spurious wave back into the domain (see Leveque, 2002,
    % p. 134 and Toro, 2009, p. 495). Also known as: zero-gradient, 
    % non-reflective, absorbing.
    methods
        %% Constructor
        function these = Transmissive(n)
            if nargin > 0
                these(1,n) = Transmissive;
            end
        end
    end
    methods(Access = protected)
        %% Enforce (scalar, override)
        function apply_scalar(this,~,~,isLeft)
            if isLeft
                this.ghostElement.states = this.boundElement.statesL;
            else
                this.ghostElement.states = this.boundElement.statesR;
            end
            this.ghostElement.interpolateStateAtEdges
        end
    end
end