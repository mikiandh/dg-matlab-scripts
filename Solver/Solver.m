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
        timeNow
        timeStop
        courantNumber
        timeDelta
        limiter
        isTimeDeltaFixed
        iterSkip
        waitForKey
        iterationCount
        wallClockTime
        exactSolution = @(t,x) nan
    end
    properties (Access = protected)
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
                validateattributes(timeSpan,{'numeric'},{'vector','nondecreasing','nonnegative'})
                this.physics = physics;
                this.timeNow = timeSpan(1);
                this.timeStop = timeSpan(end);
                % Initialize an input parser:
                p = inputParser;
                p.KeepUnmatched = true;
                % Name-value arguments:
                addParameter(p,'courantNumber',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','finite'}))
                addParameter(p,'timeDelta',[],@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','finite'}))
                addParameter(p,'limiter',Limiter,@(x)validateattributes(x,{'Limiter'},{}))
                addParameter(p,'exactSolution',@(t,x) nan,@(x)validateattributes(x,{'function_handle'},{}))
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
            addParameter(p,'limiter',this.limiter,@(x)validateattributes(x,{'Limiter'},{}))
            addParameter(p,'initialCondition',@(x) this.exactSolution(this.timeNow,x),@(x)validateattributes(x,{'function_handle'},{}))
            % Parse the inputs:
            parse(p,varargin{:})
            % Initialize solution:
            this.iterationCount = 0; tic
            for element = mesh.elements
                element.basis.(p.Results.method)(element,p.Results.initialCondition)
                element.interpolateStateAtEdges
            end
            % Update ghost elements:
            mesh.boundaries.apply(this.physics,this)
            % Limit solution:
            p.Results.limiter.applyInitial(mesh,this)
            this.limiter.takeSnapshot(mesh)
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
            this.monitor.initialize(mesh)
            if this.waitForKey
                fprintf('Paused. Press key to continue...');
                pause
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
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
        function CFL = optimizeCFL(this,space)
            % This function estimates the maximum allowable Courant number for a given
            % combination of time and space discretization methods. See 'Examples III'
            % for more details.
            %
            function [c,ceq] = nonlcon(g)
                c = max(max(g)) - 1;
                ceq = [];
            end
            CFL = 10; % Courant number seed
            Nx = 32; % number of patches (only affects resolution)
            theta = -1i*MWA_eigen_full(Nx,space.degree,space,1); % kMod = 1i*eigenvalue
            problem.options = optimoptions('fmincon','Display','notify-detailed');
            problem.solver = 'fmincon';
            problem.objective = @(CFL) - sum(sum(abs(this.amplificationFactorFun(theta*CFL)),2)); % maximize some norm of amplitude factors, all eigenmodes at once
            problem.x0 = CFL;
            problem.lb = 0;
            problem.nonlcon = @(CFL) nonlcon(abs(this.amplificationFactorFun(theta*CFL)));
            CFL = fmincon(problem);
        end
        %% Solver status
        function info = getInfo(this)
            % Returns the status of the time marching routine.
            if this.isTimeDeltaFixed
                info = sprintf('\\varsigma = %.3g, \\Deltat = %.3g (fixed)',this.courantNumber,this.timeDelta);
            else
                info = sprintf('\\varsigma = %.3g (fixed), \\Deltat = %.3g',this.courantNumber,this.timeDelta);
            end
            info = sprintf('SSP RK%d(%d); t = %.4g, %s, iter = %d, WCT = %.2f s',...
                this.order,this.stageCount,this.timeNow,info,this.iterationCount,this.wallClockTime);
        end
    end
    methods (Access = protected)
        %% Single step forward
        function STOP = stepForward(this,mesh)
            STOP = false;
            this.iterationCount = this.iterationCount + 1;
            this.timeNow = this.timeNow + this.timeDelta;
            % Check solution boundedness:
            states = cell2mat({mesh.elements.states});
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
                % Apply limiter after each stage:
                this.limiter.applyStage(mesh,this)
            end
            % Apply limiter after a full time step (one additional time):
            this.limiter.applyStep(mesh,this)
            % Plot solution:
            if STOP || ~mod(this.iterationCount,this.iterSkip)
                this.monitor.refresh(mesh)
                this.limiter.takeSnapshot(mesh)
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