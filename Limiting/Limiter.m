classdef Limiter < handle
    %
    % Base class of all DG limiters, i.e. those which can be applied to a
    % mesh containing one or more DG patches in which p > 0. Any valid 
    % limiter must inherit this class. Implements the default (empty)
    % limiter.
    %
    % Default limiter: applies its sensor, but performs no limiting.
    %
    properties (Access = protected)
        % Any limiter can be augmented with a sensor, used to exclude
        % certain elements from its action.
        sensor = Sensor;
    end
    methods
        %% Factory constructor
        function this = Limiter(varargin)
            % Given a number of input name-value pairs, decides which
            % limiter to instantiate (the default one, or one of its 
            % children) and returns it.
            %
            %%% TO DO %%%
        end
        %% Default limiting
        function apply(this,mesh,varargin)
            % Applies its sensor but does no limiting (empty limiter).
            this.sensor.apply(mesh,varargin);
        end
    end
end