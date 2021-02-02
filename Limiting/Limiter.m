classdef Limiter < handle & matlab.mixin.Heterogeneous
    %
    % Base class of all DG limiters, i.e. those which can be applied to a
    % mesh containing one or more DG patches in which p > 0. Any valid
    % limiter must inherit this class. Implements the default (empty)
    % limiter.
    %
    % Default limiter: applies its sensor, but performs no limiting.
    %
    properties (SetAccess = protected)
        % Number of times this limiter has been applied since solver
        % initialization.
        applyCount = 0;
        % Fraction of queried DOFs that have been limited since solver
        % initialization (moving average).
        cumulativeActivationRatio = 0;
        % Any limiter can be augmented with a sensor, used to exclude
        % certain elements from its action.
        sensor
        % Position of the current one with respect to a broader sequence of
        % limiters.
        priority
    end
    properties (Access = protected)
        % Cell array of instantaneous limiter activation status flags 
        % (cell: time instant, oldest to newest; column: element;
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
        %% Apply (each stage)
        function applyStage(this,mesh,solver)
            % Method that applies a sensor and subsequenly performs 
            % limiting on all elements of a given mesh for which the 
            % "isTroubled" property has been  set to true. All state vector 
            % entries which have been modified by the limiting are 
            % indicated as such in their element's "isLimited" property.
            %
            % Apply its sensor:
            this.sensor.apply(mesh,solver,this.priority);
            % Reset isLimited fields:
            for element = mesh.elements
                element.isLimited(:,:,this.priority) = false(size(element.states));
            end
        end
        %% Apply (full time-step)
        function applyStep(~,varargin)
            % Intentionally does nothing (default case).
        end
        %% Apply (initialization)
        function applyInitial(this,mesh,solver)
            % Get physics and do as every other step (default case).
            this.physics = solver.physics;
            this.applyStage(mesh,solver)
        end
        %% Name (scalar)
        function name = getName(this)
            name = sprintf('%s (%.1f%%) + %s (%.1f%%)',...
                class(this.sensor),...
                this.sensor.cumulativeActivationRatio*100,...
                class(this),...
                this.cumulativeActivationRatio*100);
            name = strrep(name,'Sensor (100.0%)','no sensor');
            name = strrep(name,'Limiter (0.0%)','no limiter');
        end
        %% Update cumulative statistics
        function updateStats(this,mesh)
            % Measure instantaneous activation ratios (sensor and limiter):
            troubledDofs = 0;
            limitedDofs = 0;
            for element = mesh.elements
                troubledDofs = troubledDofs + element.isTroubled(:,:,this.priority)*numel(element.isLimited(:,:,this.priority));
                limitedDofs = limitedDofs + sum(sum(element.isLimited(:,:,this.priority)));
            end
            troubledDofs = troubledDofs/(mesh.dofCount*this.physics.equationCount);
            limitedDofs = limitedDofs/(mesh.dofCount*this.physics.equationCount);
            % Update 'is troubled' moving average:
            this.sensor.cumulativeActivationRatio = this.sensor.cumulativeActivationRatio*this.sensor.applyCount + troubledDofs;
            this.sensor.applyCount = this.sensor.applyCount + 1;
            this.sensor.cumulativeActivationRatio = this.sensor.cumulativeActivationRatio/this.sensor.applyCount;
            % Update 'is limited' moving average:
            this.cumulativeActivationRatio = this.cumulativeActivationRatio*this.applyCount + limitedDofs;
            this.applyCount = this.applyCount + 1;
            this.cumulativeActivationRatio = this.cumulativeActivationRatio/this.applyCount;
        end
        %% Record status
        function takeSnapshot(this,mesh)
            % Stores the limiter activation status from a mesh and adds it
            % to the cell array of stored "limiter snapshots". Only stores
            % the 1st system component. Also calls its sensor's analogue.
            %
            for k = mesh.elementCount:-1:1
                aux{k,1} = sum(mesh.elements(k).isLimited(1,:,this.priority));
            end
            this.snapshots = [this.snapshots aux];
            this.sensor.takeSnapshot(mesh,this.priority);
        end
        %% Display history
        function viewTimeline(this)
            % Visualizes, in a 3D bar chart, every snapshot of this limiter
            % and its sensor at once. Quick and dirty.
            %
            figure
            subplot(2,1,1)
            b = bar3(cell2mat(this.snapshots')');
            for k = 1:length(b)
                zdata = b(k).ZData;
                b(k).CData = zdata;
                b(k).FaceColor = 'interp';
            end
            xlabel('Iteration')
            set(gca,'XTickLabel',xticks-1)
            ylabel('Cell index')
            zlabel('Number of limited modes')
            view(-90,90)
            colormap(distinguishable_colors(max(zlim),{'w','k'}))
            title('Limiter')
            subplot(2,1,2)
            this.sensor.viewTimeline
            title('Sensor')
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
    methods (Sealed)
        %% Information (heterogeneous array)
        function info = getInfo(these)
            info = these(1).getName;
            for this = these(2:end)
                info = sprintf('%s + %s',info,this.getName);
            end
        end
        %% Reset priorities (heterogeneous array)
        function resetPriorities(these)
            priorities = num2cell(1:numel(these));
            [these.priority] = priorities{:};
        end
        %% Reset physics (heterogeneous array)
        function resetPhysics(these,solver)
            [these.physics] = deal(solver.physics);
        end
        %% Reset activation statistics (heterogeneous array)
        function resetStats(these)
            for this = these
                this.sensor.applyCount = 0;
                this.sensor.cumulativeActivationRatio = 0;
                this.applyCount = 0;
                this.cumulativeActivationRatio = 0;
            end
        end
    end
end