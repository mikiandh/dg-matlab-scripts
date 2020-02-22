classdef RightBoundary < Edge
    methods
        %% Constructor
        function edge = RightBoundary(x,eL)
            if nargin > 0
                edge.coord = x;
                edge.elementL = eL;
                eL.edgeR = edge;
            end
        end
        %% Update local time-step size
        function computeTimeDeltas(this,waveSpeeds)
            this.elementL.localTimeDelta = min([this.elementL.localTimeDelta; - this.elementL.dx./waveSpeeds(waveSpeeds < 0)]); %/(2*this.elementL.basis.degree + 1)]);
        end
    end
end