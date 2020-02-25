classdef Sensor < handle
    %
    % Base class of all DG sensors, i.e. those which can be applied to a
    % mesh containing one or more DG patches. Any valid sensor must inherit
    % this class.
    %
    methods
        %% Apply (default)
        function apply(~,mesh,~)
            % Method that sets the "isTroubled" property of certain
            % elements in a given mesh to true. Any limiter will only act
            % on elements that trigger its sensor - i.e. for which this
            % method has set "isTroubled" to true.
            %
            % Mark all elements as troubled (i.e. to be limited):
            set(mesh.elements(2:end-1),'isTroubled',true); %%% FIXME: find a better way to deal with BCs %%%
            % Exclude any p == 0 cells:
            elements = findobj(mesh.elements,'dofCount',1);
            set(elements,'isTroubled',false);
        end
    end
end