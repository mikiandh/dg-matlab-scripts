classdef NoSensor < Sensor
    %
    % Default sensor: all p > 0 patches are considered troubled.
    %
    methods
        %% Default sensor
        function apply(~,mesh,varargin)
            % Mark all elements as troubled (i.e. to be limited):
            set(mesh.elements(2:end-1),'isTroubled',true); %%% FIXME: find a better way to deal with BCs %%%
            % Exclude any p == 0 cells (cannot be limited):
            elements = findobj(mesh.elements,'dofCount',1);
            set(elements,'isTroubled',false);
        end
    end
end