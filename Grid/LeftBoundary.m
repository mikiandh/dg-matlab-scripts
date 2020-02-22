classdef LeftBoundary < Edge
    methods
        %% Constructor
        function edge = LeftBoundary(x,eR)
            if nargin > 0
                edge.coord = x;
                edge.elementR = eR;
                eR.edgeL = edge;
            end
        end
        %% Update local time-step size
        function computeTimeDeltas(this,waveSpeeds)
            this.elementR.localTimeDelta = min([this.elementR.localTimeDelta; this.elementR.dx./waveSpeeds(waveSpeeds > 0)]); %/(2*this.elementR.basis.degree + 1)]);
        end
    end
end