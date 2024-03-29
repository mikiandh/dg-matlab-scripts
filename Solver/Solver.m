classdef Solver < matlab.mixin.SetGet
    % Superclass to be inherited by all explicit time-integrators of the
    % SSP RK family. Interface to various methods.
    properties (Abstract,Constant)
        order
        stageCount
        amplificationFactorFun
    end
    properties
        physics
        timeStart
        timeNow
        timeStop
        courantNumber
        timeDelta
        limiters
        isTimeDeltaFixed
        iterSkip
        waitForKey
        iterationCount
        wallClockTime
        exactSolution
    end
    properties (SetAccess = protected)
        stageNow
        initialCondition = @(x) nan
        monitor
    end
    methods (Abstract)
        %% RK stage update
        applyStage(this,element,stage)
    end
    methods
        %% Constructor
        function this = Solver(physics,timeSpan,varargin)
            if nargin > 0
                % Required inputs:
                validateattributes(physics,{'Physics'},{'scalar'})
                validateattributes(timeSpan,{'numeric'},{'vector','nondecreasing','nonnegative','numel',2})
                this.physics = physics;
                this.timeStart = timeSpan(1);
                this.timeStop = timeSpan(end);
                this.timeNow = this.timeStart;
                % Initialize an input parser:
                p = inputParser;
                p.KeepUnmatched = true;
                % Name-value arguments:
                addParameter(p,'courantNumber',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','finite'}))
                addParameter(p,'timeDelta',[],@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','finite'}))
                addParameter(p,'limiters',Limiter,@(x)validateattributes(x,{'Limiter'},{}))
                addParameter(p,'exactSolution',@(t,x) nan(physics.equationCount,numel(x)),@(x)validateattributes(x,{'function_handle'},{}))
                addParameter(p,'iterSkip',0,@(x)validateattributes(x,{'numeric'},{'integer'}))
                addParameter(p,'waitForKey',false,@(x)validateattributes(x,{'logical'},{'scalar'}))
                % Parse the inputs:
                parse(p,varargin{:})
                set(this,fieldnames(p.Results)',struct2cell(p.Results)')
                % If timeDelta has been given a value, keep it fixed:
                if any(string(p.UsingDefaults) == "timeDelta")
                    this.isTimeDeltaFixed = false;
                else
                    this.isTimeDeltaFixed = true;
                end
                % Initialize a monitor for this solver instance:
                this.monitor = Monitor(this,varargin{:});
            end
        end
        %% Initialize solution
        function initialize(this,mesh,varargin)
            % Initializes a given mesh by calling a given projection type
            % (L2 projection by default) on each element. Exact solution,
            % limiter and sensor can be overriden (those of the solver are
            % used by default).
            %
            % Initialize an input parser:
            p = inputParser;
            % Optional arguments:
            addParameter(p,'method','project',@(x)validateattributes(x,{'char'},{}))
            addParameter(p,'limiters',this.limiters,@(x)validateattributes(x,{'Limiter'},{}))
            addParameter(p,'initialCondition',@(x) this.exactSolution(this.timeNow,x),@(x)validateattributes(x,{'function_handle'},{}))
            % Parse the inputs:
            parse(p,varargin{:})
            % Inform custom limiters of their order:
            p.Results.limiters.resetPriorities;
            % Initialize solution:
            this.iterationCount = 0; tic
            for element = mesh.elements
                element.basis.(p.Results.method)(element,p.Results.initialCondition)
            end
            mesh.elements.interpolateStateAtEdges
            % Update ghost elements:
            mesh.boundaries.apply(this.physics,this)
            % Limit initial solution (possibly using custom limiters):
            for limiter = p.Results.limiters
                limiter.applyInitial(mesh,this)
                limiter.takeSnapshot(mesh)
            end
            p.Results.limiters.resetStats
            p.Results.limiters.updateStats(mesh)
            this.wallClockTime = toc;
            % Initialize residuals:
            mesh.computeResiduals(this.physics,this)
            % Initialize time-step size:
            if this.isTimeDeltaFixed
                this.updateCourantNumber(mesh)
            else
                this.updateTimeDelta(mesh)
            end
            % Initialize solver monitor (shows initial condition):
            this.monitor.initialize(mesh,p.Results.initialCondition)
            if this.waitForKey
                fprintf('Paused. Press key to continue...');
                pause
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
            end
            % Adjust limiter stuff, if necessary.
            if any(strcmpi(p.UsingDefaults,'limiters'))
                % Default case, no extra adjustments needed.
            else
                % Inform each owned limiter of its position in the sequence:
                this.limiters.resetPriorities
                % And also of the current physics:
                this.limiters.resetPhysics(this)
                % And also reset its statistics, just in case:
                this.limiters.resetStats
                % Reset and preallocate activation flags (in case limiters were not default):
                for element = mesh.elements
                    element.isTroubled = repmat(element.isTroubled(:,:,1),1,1,numel(this.limiters));
                    element.isLimited = repmat(element.isLimited(:,:,1),1,1,numel(this.limiters));
                end
            end
        end
        %% Advance solution in time
        function launch(this,mesh)
            % Start and drive the time marching routine, from timeNow
            % until timeStop of this solver instance.
            %
            tic
            STOP = this.timeNow == this.timeStop;
            if this.isTimeDeltaFixed
                while ~STOP
                    STOP = this.stepForward(mesh);
                    this.updateCourantNumber(mesh)
                    this.wallClockTime = toc;
                end
            else
                while ~STOP
                    STOP = this.stepForward(mesh);
                    this.updateTimeDelta(mesh)
                    this.wallClockTime = toc;
                end
            end
            if ~STOP
                warning('Solver stopped prematurely.')
            end
        end
        %% Estimate maximum stable Courant number
        function CFL = optimizeCFL(this,space,varargin)
            % This function estimates the maximum allowable Courant number for a given
            % combination of time and space discretization methods. See 'Examples III'
            % for more details.
            %
            function [c,ceq] = nonlcon(x)
                c = max(abs(this.amplificationFactorFun(theta*x))) - 1;
                ceq = [];
            end
            theta = reshape(space.getFourierFootprint(varargin{:}),1,[]); % 1D array
            CFL = fmincon(@(CFL) -CFL, 1,[],[],[],[],0,[],@nonlcon,optimoptions('fmincon','Display','notify'));
        end
        %% Solver status
        function info = getInfo(this)
            % Returns the status of the time marching routine.
            if this.isTimeDeltaFixed
                info = sprintf('\\varsigma = %.3g, \\Deltat = %.3g (fixed)',this.courantNumber,this.timeDelta);
            else
                info = sprintf('\\varsigma = %.3g (fixed), \\Deltat = %.3g',this.courantNumber,this.timeDelta);
            end
            info = sprintf('SSP RK%d(%d); t = %.4g \\in [%.4g,%.4g], %s, iter = %d, WCT = %.2f s',...
                this.order,this.stageCount,this.timeNow,this.timeStart,this.timeStop,info,this.iterationCount,this.wallClockTime);
        end
        %% Export solution
        function writeSolutionToFile(this,fileNameRoot,n)
            % Exports the solution (numerical and exact; also initial
            % condition) into a .dat file.
            %
            % Arguments
            %  fileNameRoot: root of the name the file will have
            %  n: write only one of every n samples (4 by default)
            %
            % Setup file:
            fileName = sprintf('%s.dat',fileNameRoot);
            fileID = fopen(fileName,'wt');
            % Pass on request to the monitor:
            if nargin < 3
                n = 4;
            end
            try
                this.monitor.writeInfo(fileID)
                this.monitor.writeSolution(fileID,n)
            catch me
                warning(getReport(me))
            end
            % Return file to a safe state:
            fclose(fileID);
        end
        %% Export limiter & sensor
        function writeLimiterToFile(this,fileNameRoot)
            % Exports the current limiter and sensor activation status into
            % a .dat file.
            %
            % Arguments
            %  fileNameRoot: root of the name the file will have
            %
            % Setup file:
            fileName = sprintf('%s.dat',fileNameRoot);
            fileID = fopen(fileName,'wt');
            % Pass on request to the monitor:
            try
                this.monitor.writeInfo(fileID)
                this.monitor.writeLimiter(fileID)
            catch me
                warning(getReport(me))
            end
            % Return file to a safe state:
            fclose(fileID);
        end
        %% Export solution points
        function writePointsToFile(this,fileNameRoot,ptName)
            % Exports samples of the discrete solution at selected points.
            %
            % Arguments
            %  fileNameRoot: root of the name the file will have
            %  ptName: which type of points to export (will be appended to
            %          the file name)
            %
            % Setup file:
            fileName = sprintf('%s_%s.dat',fileNameRoot,ptName);
            fileID = fopen(fileName,'wt');
            % Pass on request to the monitor:
            try
                this.monitor.writeInfo(fileID)
                this.monitor.writePoints(fileID,ptName)
            catch me
                warning(getReport(me))
            end
            % Return file to a safe state:
            fclose(fileID);
        end
    end
    methods (Access = protected)
        %% Single step forward
        function STOP = stepForward(this,mesh)
            STOP = false;
            this.iterationCount = this.iterationCount + 1;
            this.timeNow = this.timeNow + this.timeDelta;
            % Check solution boundedness:
            states = [mesh.elements.states];
            if any(isinf(states(:))) || any(isnan(states(:)))
                warning('Inf or NaN state(s) detected at t = %.4f.',this.timeNow)
                STOP = true;
                return
            end
            % Check stop criterion:
            if this.timeNow >= this.timeStop
                this.timeDelta = this.timeStop - this.timeNow + this.timeDelta;
                this.timeNow = this.timeStop;
                STOP = true;
            end
            % Advance one time-step:
            this.stageNow = 0; % reset stage counter
            while this.stageNow < this.stageCount
                % Advance stage counter:
                this.stageNow = this.stageNow + 1;
                % Evaluate solution residuals:
                mesh.computeResiduals(this.physics,this)
                % Advance solution by one stage:
                for element = mesh.elements
                    this.applyStage(element)
                end
                % Apply limiter(s) after each stage:
                for limiter = this.limiters
                    limiter.applyStage(mesh,this)
                end
                this.limiters.updateStats(mesh)
            end
            % Apply limiter(s) after a full time step (one additional time):
            for limiter = this.limiters
                limiter.applyStep(mesh,this)
            end
            % Plot solution:
            if STOP || ~mod(this.iterationCount,this.iterSkip)
                this.monitor.refresh(mesh)
                for limiter = this.limiters
                    limiter.takeSnapshot(mesh)
                end
                if this.waitForKey
                    fprintf('Paused. Press key to continue...');
                    pause
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
                end
            end
        end
        %% Update Courant number (fixed timeDelta)
        function updateCourantNumber(this,mesh)
            % Set global Courant number to smallest local one.
            %
            this.courantNumber = inf;
            for element = mesh.elements
                this.courantNumber = min([this.courantNumber element.localTimeDelta]);
            end
            this.courantNumber = this.timeDelta/this.courantNumber;
        end
        %% Update time-step size (fixed Courant number)
        function updateTimeDelta(this,mesh)
            % Set global time-step size to smallest local one.
            %
            this.timeDelta = inf;
            for element = mesh.elements
                this.timeDelta = min([this.timeDelta element.localTimeDelta]);
            end
            this.timeDelta = this.courantNumber*this.timeDelta;
        end
    end
end