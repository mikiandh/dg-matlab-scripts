classdef NoLimiter < Limiter
    %
    % Default limiter: applies its sensor, but performs no limiting. Sets
    % all isLimited elements to false (i.e. no limiting was performed).
    %
    methods
        %% Default limiting
        function apply(this,mesh,varargin)
            % Apply its sensor:
            this.sensor.apply(mesh,varargin);
            % Reset isLimited fields:
            for element = mesh.elements
                element.isLimited = false(size(element.states));
            end
        end
    end
end