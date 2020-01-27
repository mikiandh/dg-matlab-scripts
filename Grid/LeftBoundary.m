classdef LeftBoundary < Edge
    methods
        %% Constructor
        function edge = LeftBoundary(x,eR)
            if nargin > 0
                edge.coord = x;
                edge.elementR = eR;
                %edge.leftNodeRelativePosition = -inf;
                %edge.rightNodeRelativePosition = ...
                %    eR.mapFromReference(eR.DG.nodeCoords(1)) - edge.coord;
            end
        end
        %% Update local time-step size
        function computeTimeDeltas(this,waveSpeeds)
            this.elementR.localTimeDelta = min([this.elementR.localTimeDelta; this.elementR.dx./waveSpeeds(waveSpeeds > 0)]); %/(2*this.elementR.basis.degree + 1)]);
        end
    end
end