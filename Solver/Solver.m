classdef Solver < handle
    properties (Abstract,Constant)
        order
        stageCount
        amplificationFactorFun
    end
    properties
        family = 'SSP RK'
        courantNumber
        iterationCount
        timeNow
        timeDelta
        isTimeDeltaFixed
        timeStop
        physics
        limiter
        plotData
    end
    methods (Abstract)
        %% RK stage update
        applyStage(this,element,stage)
    end
    methods
        %% Constructor
        function this = Solver(t,tEnd,CFL,dt,limiter)
            if nargin > 0
                this.timeNow = t;
                this.timeStop = tEnd;
                this.courantNumber = CFL;
                this.timeDelta = dt;
                if isempty(limiter)
                    this.limiter = Limiter;
                else
                    this.limiter = limiter;
                end
            end
        end
        %% Drive the integration
        function STOP = launch(this,mesh,replotIters,solutionFun)
            % Safety check:
            if isempty(this.timeDelta)
                this.isTimeDeltaFixed = false;
                this.timeDelta = 0;
            elseif ~isempty(this.timeDelta)
                this.isTimeDeltaFixed = true;
            else
                error('Specify time-step size and/or CFL.')
            end
            % Initialize some stuff:
            STOP = false;
            this.iterationCount = 0;
            this.initializePlot(mesh,solutionFun);
            this.refreshPlot(mesh);
            % Time marching routine:
            if this.isTimeDeltaFixed
                while ~STOP
                    STOP = this.stepForward(mesh,replotIters);
                    this.updateCourantNumber(mesh);
                end
            else
                while ~STOP
                    STOP = this.stepForward(mesh,replotIters);
                    this.updateTimeDelta(mesh);
                end
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
    end
    methods (Access = protected)
        %% Single step forward
        function STOP = stepForward(this,mesh,replotIters)
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
            for stage = 1:this.stageCount
                % Evaluate solution residuals:
                mesh.computeResiduals;
                % Advance solution by one stage:
                for element = mesh.elements
                    this.applyStage(element,stage);
                end
                % Apply limiter:
                this.limiter.apply(mesh,this);
            end
            % Plot solution:
            if STOP || ~mod(this.iterationCount,replotIters)
                this.refreshPlot(mesh);
                this.limiter.takeSnapshot(mesh); % could be expensive
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
        %% Initialize plot data
        function initializePlot(this,mesh,fun)
            figDims = [800 100 700, min(mesh.physics.equationCount*400,800)];
            this.plotData.fig = figure('Renderer','painters','Position',figDims);
            this.plotData.exactFun = fun;
            this.plotData.cmap = distinguishable_colors(mesh.elementCount,{'w','k'});
            this.plotData.xGlobal = linspace(mesh.elements(1).xL,mesh.elements(end).xR,1000);
            this.plotData.ylims = repmat([0 1.2],mesh.physics.equationCount,1);
            disc = mesh.bases(1);
            this.plotData.space = strrep(class(disc),'_','-');
            this.plotData.limiter = this.limiter.getInfo;
            % Deduce method name:
            if isa(disc,'FR')
                aux = sprintf('%s(%s)',this.plotData.space,num2str(disc.param));
            elseif disc.isHybrid
                if isempty(disc.smoothness)
                    aux = sprintf('%s, N_\\Sigma = %d',this.plotData.space,disc.nonzeroSpanCount);
                else
                    if isnan(disc.smoothness) || isinf(disc.smoothness)
                        varkappa = disc.degree - 1;
                    else
                        varkappa = disc.smoothness;
                    end
                    aux = sprintf('%s, N_\\Sigma = %d, C^{%d}',this.plotData.space,disc.nonzeroSpanCount,varkappa);
                end
            else
                aux = sprintf('%s',this.plotData.space);
            end
            % Append limiter info:
            aux = sprintf('%s; %s',aux,this.plotData.limiter);
            % Set up plot title:
            this.plotData.space =...
                sprintf('%s;  p \\in [%d,%d], N_\\Omega = %d, N = %d',...
                aux,mesh.minDegree,mesh.maxDegree,mesh.elementCount,mesh.dofCount);
            this.plotData.time = sprintf('%s%d(%d)',this.family,this.order,this.stageCount);
        end
        %% Plot all state vector components
        function refreshPlot(this,mesh)
            if this.isTimeDeltaFixed
                aux = sprintf('\\varsigma = %.3g, \\Deltat = %.3g (fixed)',this.courantNumber,this.timeDelta);
            else
                aux = sprintf('\\varsigma = %.3g (fixed), \\Deltat = %.3g',this.courantNumber,this.timeDelta);
            end
            aux = sprintf('; t = %.4g, %s, iter = %d',this.timeNow,aux,this.iterationCount);
            figure(this.plotData.fig)
            yGlobal = this.plotData.exactFun(this.timeNow,this.plotData.xGlobal);
            for i = 1:mesh.physics.equationCount
                subplot(mesh.physics.equationCount,1,i)
                plot(this.plotData.xGlobal,yGlobal(i,:),'-.k')
                hold on
                k = 1;
                for element = mesh.elements
                    samples = 1+(2*element.basis.basisCount)^2;
                    %x = sort([linspace(element.xL,element.xR,samples) element.getNodeCoords']);
                    x = sort(linspace(element.xL,element.xR,samples));
                    y = element.interpolateStateAtCoords(x);
                    plot(x,y(i,:),'-','Color',this.plotData.cmap(k,:))
                    if element.basis.isNodal
                        plot(element.getNodeCoords,element.states(i,:),...
                            '.','Color',this.plotData.cmap(k,:),'MarkerSize',10)
                    end
                    if element.basis.isHybrid
                        x = element.getBreakCoords;
                        y = element.interpolateStateAtCoords(x);
                        hLine = scatter(x,y(i,:),'LineWidth',0.25);
                        set(hLine.MarkerHandle,'Style','vbar');
                        set(hLine,'MarkerFaceColor',this.plotData.cmap(k,:));
                        set(hLine,'MarkerEdgeColor',this.plotData.cmap(k,:));
                    end
                    k = k + 1;
                end
                hold off
                xlabel('x')
                ylabel(['q_' num2str(i) '(t,x)'])
                setFancyPlotSettings3
                yl = ylim;
                this.plotData.ylims(i,:) = [min(yl(1),this.plotData.ylims(i,1))...
                    max(yl(2),this.plotData.ylims(i,2))];
                ylim(this.plotData.ylims(i,:))
                if i == 1
                    title({...
                        class(mesh.physics),...
                        this.plotData.space,...
                        strcat(this.plotData.time,aux)...
                        });
                end
            end
            drawnow limitrate
        end
    end
end