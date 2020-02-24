classdef Sensor < handle
    %
    % Base class of all DG sensors, i.e. those which can be applied to a
    % mesh containing one or more DG patches. Any valid sensor must inherit
    % this class.
    %
    % Default sensor: all p > 0 patches are sensed as possibly troubled.
    %
    methods
        %% Factory constructor
        function this = Sensor(varargin)
            % Given a number of input name-value pairs, decides which
            % sensor to instantiate and returns it.
            %
            %%% Default sensor (empty) %%%
        end
        %% Default sensing
        function apply(this,mesh,varargin) %#ok<INUSL>
            % Turns off limiting at all p == 0 cells:
            elements = findobj(mesh.elements,'dofCount',1);
            set(elements,'isSensed',false);
        end
    end
end