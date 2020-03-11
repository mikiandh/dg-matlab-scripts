classdef Limiter < handle
    %
    % Base class of all DG limiters, i.e. those which can be applied to a
    % mesh containing one or more DG patches in which p > 0. Any valid
    % limiter must inherit this class. Implements the default (empty)
    % limiter.
    %
    % Default limiter: applies its sensor, but performs no limiting.
    %
    properties (SetAccess = protected)
        % Any limiter can be augmented with a sensor, used to exclude
        % certain elements from its action.
        sensor
    end
    properties (Access = protected)
        % The activation status of the limiter is accumulated in this cell
        % array (cell: time instant, oldest to newest; column: element;
        % row: equation component).
        snapshots
        % In order to project to/from characteristic variables, any limiter
        % will need to be aware of the physics being solved.
        physics
    end
    methods
        %% Constructor (default)
        function this = Limiter(varargin)
            % Instantiates the default limiter (which does no limiting)
            % with a sensor parsed from a given list of name-value pairs.
            %
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'Sensor',Sensor,@(x) validateattributes(x,{'Sensor'},{}));
            % Parse the input sensor:
            parse(p,varargin{:});
            this.sensor = p.Results.Sensor;
        end
        %% Apply (default)
        function apply(this,mesh,varargin)
            % Method that applies a sensor and subsequenly performs 
            % limiting on all elements of a given mesh for which the 
            % "isTroubled" property has been  set to true. All state vector 
            % entries which have been modified by the limiting are 
            % indicated as such in their element's "isLimited" property.
            %
            % Retrieve physics:
            this.physics = mesh.physics;
            % Apply its sensor:
            this.sensor.apply(mesh,varargin{:});
            % Reset isLimited fields:
            for element = mesh.elements
                element.isLimited = false(size(element.states));
            end
        end
        %% Information
        function info = getInfo(~)
            info = 'no limiting';
        end
        %% Record status
        function takeSnapshot(this,mesh)
            % Stores the limiter activation status from a mesh and adds it
            % to the cell array of stored "limiter snapshots". Only stores
            % the 1st system component.
            %
            aux = {mesh.elements.isLimited};
            aux = cellfun(@(x) sum(x(1,:),2),aux,'UniformOutput',false);
            this.snapshots = [this.snapshots aux'];
        end
        %% Display history
        function viewTimeline(this)
            % Visualizes, in a 3D bar chart, every snapshot of the limiter
            % at once. Quick and dirty.
            %
            b = bar3(cell2mat(this.snapshots')');
            colorbar
            for k = 1:length(b)
                zdata = b(k).ZData;
                b(k).CData = zdata;
                b(k).FaceColor = 'interp';
            end
            xlabel('Time level')
            ylabel('Cell index')
            zlabel('Number of limited modes')
            view(-90,90)
        end
    end
    methods (Static)
        %% Singleton constructor
        function limiter = get(name,varargin)
            % Returns a limiter instance of the given class name, with the
            % passed list of name-value argument pairs.
            %
            if isempty(name)
                name = 'none';
            end
            switch name
                case {'Limiter','none',''}
                    limiter = Limiter(varargin{:});
                case {'TVB','TVD','Shu'}
                    limiter = TVB(varargin{:});
                case {'BDF','Biswas'}
                    limiter = BDF(varargin{:});
                case {'BSB','Burbeau'}
                    limiter = BSB(varargin{:});
                case {'Krivodonova','Kriv'}
                    limiter = Krivodonova(varargin{:});
                otherwise
                    error('Limiter name unknown.')
            end
        end
    end
end