classdef Sensor < handle
    %
    % Base class of all DG sensors, i.e. those which can be applied to a
    % mesh containing one or more DG patches. Any valid sensor must inherit
    % this class.
    %
    properties (Access = protected)
        % The activation status of the sensor is accumulated in this cell
        % array (cell: time instant, oldest to newest; column: element)
        snapshots
    end
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
        %% Record status
        function takeSnapshot(this,mesh)
            % Stores the sensor activation status from a mesh and adds it
            % to the cell array of stored "sensor snapshots".
            %
            this.snapshots = [this.snapshots {[mesh.elements.isTroubled]}];
        end
        %% Display history
        function viewTimeline(this)
            % Visualizes, in a 3D bar chart, every snapshot of this sensor
            % at once. Quick and dirty.
            %
            b = bar3(cell2mat(this.snapshots')');
            for k = 1:length(b)
                zdata = b(k).ZData;
                b(k).CData = zdata;
                b(k).FaceColor = 'interp';
            end
            xlabel('Iteration')
            set(gca,'XTickLabel',xticks-1)
            ylabel('Cell index')
            zlabel('Cell is troubled')
            view(-90,90)
        end
    end
end