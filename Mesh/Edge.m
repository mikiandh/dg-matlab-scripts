classdef Edge < matlab.mixin.SetGet
    properties
        coord
        elementL
        elementR
    end
    methods
        %% Constructor (vector)
        function these = Edge(elements)
            % Instantiates an array of (internal) edges. Sets connectivity 
            % between neighbors (edge-element and element-element).
            if nargin > 0
                these(1,numel(elements)-1) = Edge;
                [these(1:end-1).coord] = elements(2:end-1).xL;
                [these(2:end).coord] = elements(2:end-1).xR;
                % Edge-to-elements connectivity:
                aux = num2cell(elements);
                [these.elementL] = aux{1:end-1};
                [these.elementR] = aux{2:end};
                % Element-to-edges connectivity:
                edges = num2cell(these);
                [elements(2:end).edgeL] = edges{:};
                [elements(1:end-1).edgeR] = edges{:};
                % Element-to-elements connectivity:
                [elements(2:end).elementL] = aux{1:end-1};
                [elements(1:end-1).elementR] = aux{2:end};
            end
        end
        %% Riemann flux at edges
        function computeFlux(these,physics)
            % Computes and sets numerical fluxes (i.e. solutions of the
            % Riemann problem at an edge) on both elements sharing each of 
            % these edges, including a conventional sign. Also sets a local
            % tim-step size based on each Riemann problem's characteristic 
            % speeds. Vector input.
            for this = these
                [flux,waveSpeeds] = physics.riemannFlux(...
                    this.elementL.stateR,this.elementR.stateL);
                this.elementL.riemannR = flux; % 'normal vector' = +1
                this.elementR.riemannL = -flux; % 'normal vector' = -1
                this.computeTimeDeltas(waveSpeeds);
            end
        end
    end
    methods (Access = protected)
        %% Update local time-step sizes
        function computeTimeDeltas(this,waveSpeeds)
            % Sets local time-step sizes based on characteristic speeds and
            % element sizes, such that each local Courant number is unity.
            % Does not override more conservative existing values.
            this.elementL.localTimeDelta = min([this.elementL.localTimeDelta; - this.elementL.dx./waveSpeeds(waveSpeeds < 0)]);
            this.elementR.localTimeDelta = min([this.elementR.localTimeDelta; this.elementR.dx./waveSpeeds(waveSpeeds > 0)]);
        end
    end
end