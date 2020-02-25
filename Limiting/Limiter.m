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
        sensor = NoSensor; 
    end
    properties (SetAccess = protected)
        % The activation status of the limiter is accumulated in this cell
        % array (cell: time instant, oldest to newest; column: element; 
        % row: equation component).
        %
        % To visualize, try:
        %
        %  b = bar3(cell2mat(limiter.snapshots')');
        %  colorbar
        %  for k = 1:length(b)
        %      zdata = b(k).ZData;
        %      b(k).CData = zdata;
        %      b(k).FaceColor = 'interp';
        %  end
        %
        snapshots;
    end
    methods
        %% Factory constructor
        function this = Limiter(varargin)
            % Given a list of input name-value pairs, decides which
            % limiter to instantiate and returns it.
            %
            if nargin
               
            end
        end
        %% Accumulate activation history
        function takeSnapshot(this,mesh)
            % Stores the limiter activation status from a mesh and adds it
            % to the cell array of stored "limiter snapshots".
            aux = {mesh.elements.isLimited};
            aux = cellfun(@(x) sum(x,2), aux);
            this.snapshots = [this.snapshots {aux}];
        end
    end
    methods (Abstract)
        % Method that applies a sensor and subsequenly performs limiting on
        % all elements of a given mesh for which the "isSensed" property 
        % has been  set to true. All state vector entries which have been
        % modified by the limiting are indicated as such in their element's
        % "isLimited" logical 2D array.
        apply(this,mesh,solver)
    end
end