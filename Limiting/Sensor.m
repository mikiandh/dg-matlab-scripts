classdef Sensor < handle
    %
    % Base class of all DG sensors, i.e. those which can be applied to a
    % mesh containing one or more DG patches. Any valid sensor must inherit
    % this class.
    %
    methods
        %% Factory constructor
        function this = Sensor(varargin)
            % Given a list of input name-value pairs, decides which
            % sensor to instantiate and returns it.
            %
            if nargin
               
            end
        end
    end
    methods (Abstract)
        % Method that sets the "isTroubled" property of certain elements in a
        % given mesh to true. Any limiter will only act on elements that
        % trigger its sensor - i.e. for which this method has set 
        % "isTroubled" to true.
        apply(this,mesh,solver)
    end
end