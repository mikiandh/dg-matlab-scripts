classdef Edge < handle
    properties
        coord
        elementL
        elementR
    end
    methods
        %% Constructor
        function edge = Edge(x,eL,eR)
            if nargin > 0
                edge.coord = x;
                edge.elementL = eL;
                edge.elementR = eR;
            end
        end
        %% Upwind flux at interior edge
        function computeFlux(this,physics)
            [flux,waveSpeeds] = physics.riemannFlux(...
                this.elementL.stateR,this.elementR.stateL);
            this.elementL.riemannR = flux; % 'normal vector' = +1
            this.elementR.riemannL = -flux; % 'normal vector' = -1
            this.computeTimeDeltas(waveSpeeds);
        end
        %% Update local time-step sizes
        function computeTimeDeltas(this,waveSpeeds)
            this.elementL.localTimeDelta = min([this.elementL.localTimeDelta; - this.elementL.dx./waveSpeeds(waveSpeeds < 0)]); %/(2*this.elementL.basis.degree + 1)]);
            this.elementR.localTimeDelta = min([this.elementR.localTimeDelta; this.elementR.dx./waveSpeeds(waveSpeeds > 0)]); %/(2*this.elementR.basis.degree + 1)]);
        end
    end
end