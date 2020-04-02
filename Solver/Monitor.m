classdef Monitor < handle
    % Class that handles monitoring tasks related to a particular solver
    % instance. These include plotting of the solution and computation of
    % solution and/or error norms.
    properties (Constant)
        colorSkip = {'w','k'}
        nPerDof = 6
    end
    properties (SetAccess = immutable)
        solver
        norms
        showSensor
        showLimiter
        rows
        cols
    end
    properties (Access = protected)
        % Axis limits (seed) for the solution (i.e. left column) subplots:
        ylims = [0 1] % columns: min, max (#rows is set automatically)
        % Colormaps:
        cDiscrete
        cLimiters
        cNorms
        % Handles to graphics instances:
        hTitle % handle to the title
        hFigure % handle to the figure
        hAxes % 2D array of axes (i.e. subplot) handles
        hDiscrete % 2D array of patch-wide approximate solutions (row: equation; column: cell)
        hExact % idem, exact solutions
        hNodes % idem, nodal coordinates and values
        hBreaks % idem, breakpoints
        hSensors % column array of sensor activation bar charts
        hLimiters % idem, limiter activation bar chars (stacked)
        hNorms % 2D array of instantaneous norm samples (row: equation; column: norm type)
    end
    methods
        %% Constructor
        function this = Monitor(solver,varargin)
            % Initialize parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'solver',@(x)validateattributes(x,{'Solver'},{}))
            addParameter(p,'norms',[],@(x)validateattributes(x,{'Norm'},{'vector'}))
            addParameter(p,'showSensor',false,@(x)validateattributes(x,{'logical'},{'scalar'}))
            addParameter(p,'showLimiter',false,@(x)validateattributes(x,{'logical'},{'scalar'}))
            addParameter(p,'equations',1:solver.physics.equationCount,@(x)validateattributes(x,{'numeric'},{'vector','integer','>',0,'<=',solver.physics.equationCount}))
            p.parse(solver,varargin{:});
            % Set properties from parser:
            this.solver = p.Results.solver;
            this.norms = unique(p.Results.norms);
            this.showSensor = p.Results.showSensor;
            this.showLimiter = p.Results.showLimiter;
            % Subplot indices:
            this.rows = p.Results.equations;
            this.cols = 1:1 + ~isempty(this.norms);
            % Repmat initial axis limits:
            this.ylims = repmat(this.ylims,numel(this.rows),1);
        end
        %% Initialize
        function initialize(this,mesh)
            % Initializes a figure with a 2D arrangement of subplots.
            %
            % The left column shows one solution component per row, over
            % the entire domain, at a particular instant in time. Each
            % patch (i.e. DG element, a.k.a macro-element) has a different
            % colour. Breakpoints are indicated by vertical bar markers.
            % Nodes (if any) are shown as dot markers.
            %
            % The right column shows samples over time of component-wise 
            % solution and/or error norms, over the entire domain. Each
            % norm has a different colour. Sample points are indicated by
            % cross markers.
            %
            % Initialize a figure:
            this.hFigure = figure('Renderer','painters','Position',[400 100 700*this.cols(end), min(400*this.rows(end),800)]);
            % Initialize a "supertitle" annotation:
            this.hTitle = annotation(this.hFigure,'textbox',[0 0.9 1 0.1],...
                'String',{[this.solver.physics.getInfo '; ' mesh.boundaries.getInfo],[mesh.getInfo '; ' this.solver.limiter.getInfo],this.solver.getInfo},...
                'FontWeight','bold','EdgeColor','none',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            % Preallocate graphics objects:
            this.hDiscrete = gobjects(this.rows(end),mesh.elementCount);
            this.hExact = this.hDiscrete;
            this.hNodes = this.hDiscrete;
            this.hBreaks = this.hDiscrete;
            this.hSensors = gobjects(this.rows(end),1);
            this.hLimiters = gobjects(this.rows(end),mesh.maxBasisCount-1);
            this.hNorms = gobjects(this.rows(end),numel(this.norms));
            % Initialize 2D lattice of subplot axes:
            k = 0;
            for i = this.rows
                for j = this.cols
                    k = k + 1;
                    this.hAxes(i,j) = subplot(this.rows(end),this.cols(end),k);
                    % Push downwards (make room for the title):
                    if numel(this.rows) == 1
                        set(this.hAxes(i,j),'OuterPosition',get(this.hAxes(i,j),'OuterPosition').*[1 1 1 .85])
                    elseif numel(this.rows) == 2
                        set(this.hAxes(i,j),'OuterPosition',get(this.hAxes(i,j),'OuterPosition').*[1 1 1 .90])
                    elseif numel(this.rows) == 3
                        set(this.hAxes(i,j),'OuterPosition',get(this.hAxes(i,j),'OuterPosition').*[1 1 1 .95])
                    end
                    % Add labels:
                    if j == 1
                        xlabel('x')
                        ylabel(sprintf('q_%d(t,x)',i))
                        hold on
                        setFancyPlotSettings3
                    else
                        xlabel('t')
                        ylabel(sprintf('||q_%d(t,x)||',i))
                        setFancyPlotSettings2
                    end
                end
            end
            % Initialize colormaps:
            this.cLimiters = distinguishable_colors(mesh.maxBasisCount-1,this.colorSkip);
            %%%this.cLimiters = jet(mesh.maxBasisCount-1);
            this.cDiscrete = distinguishable_colors(mesh.elementCount,this.colorSkip);
            this.cNorms = distinguishable_colors(numel(this.norms),this.colorSkip);
            % Initialize sensor bar charts (if requested):
            if this.showSensor
                for i = this.rows
                    this.hSensors(i) = bar(this.hAxes(i,1),[mesh.elements.x],[mesh.elements.isTroubled],1);
                end
                set(this.hSensors,'FaceAlpha',.05,'EdgeAlpha',.1,...
                    'EdgeColor',this.cLimiters(1,:),'FaceColor',this.cLimiters(1,:))
            end
            % Initialize limiter bar charts (if requested):
            if this.showLimiter 
                for i = this.rows
                    this.hLimiters(i,:) = bar(this.hAxes(i,1),[mesh.elements.x],this.getLimiterData(mesh,i),1,'stacked'); % fraction of 1
                    set(this.hLimiters(i,:),{'EdgeColor'},num2cell(this.cLimiters,2),{'FaceColor'},num2cell(this.cLimiters,2));
                end
                set(this.hLimiters,'FaceAlpha',.1,'EdgeAlpha',.2)
            end
            % Initialize solution vs. space plots:
            for k = 1:mesh.elementCount
                % Sample coordinates:
                x1 = mesh.elements(k).getNodeCoords';
                x2 = mesh.elements(k).getBreakCoords;
                x3 = linspace(mesh.elements(k).xL,mesh.elements(k).xR,mesh.elements(k).dofCount*this.nPerDof);
                x = unique([x3 x2 x1],'sorted');
                % Indices of sample subsets:
                j = cellfun(@(b) arrayfun(@(a) find(x == a,1),b),{x1 x2 x3},'UniformOutput',false);
                % Solution samples:
                y = mesh.elements(k).interpolateStateAtCoords(x);
                z = this.solver.exactSolution(this.solver.timeNow,x);
                for i = this.rows
                    this.hExact(i,k) = plot(this.hAxes(i,1),x,z(i,:),'LineStyle','--','Color','k');
                    this.hDiscrete(i,k) = plot(this.hAxes(i,1),x,y(i,:),'LineStyle','-','Color',this.cDiscrete(k,:));
                    this.hBreaks(i,k) = plot(this.hAxes(i,1),x2,y(i,j{2}),'LineStyle','none','Marker','+','Color',this.cDiscrete(k,:));
                    this.hNodes(i,k) = plot(this.hAxes(i,1),[nan x1],[nan y(i,j{1})],'LineStyle','none','Marker','.','Color',this.cDiscrete(k,:));
                end
            end
            % Initialize norm vs. time plots (if any):
            if ~isempty(this.norms)
                this.norms.compute(mesh,@(x) this.solver.exactSolution(this.solver.timeNow,x));
                for j = 1:numel(this.norms)
                    for i = this.rows
                        this.hNorms(i,j) = animatedline(this.solver.timeNow,this.norms(j).vals(i),'Parent',this.hAxes(i,2),...
                            'LineStyle','-','Marker','x','Color',this.cNorms(j,:),'DisplayName',char(this.norms(j)));
                    end
                end
                % Show legend:
                legend(this.hAxes(end),'Location','NorthEast')
            end
            % Redraw:
            set(this.hAxes(:,1),'XLim',[mesh.edges(1).coord mesh.edges(end).coord])
            if ~isscalar(this.cols)
                set(this.hAxes(:,2),'XLim',[this.solver.timeStart this.solver.timeStop])
            end
            %%%this.updateYLims
            drawnow
        end
        %% Refresh
        function refresh(this,mesh)
            % Replots the monitor figure with current solution and norms.
            %
            % Refresh the title:
            this.hTitle.String = {[this.solver.physics.getInfo '; ' mesh.boundaries.getInfo], [mesh.getInfo '; ' this.solver.limiter.getInfo],this.solver.getInfo};
            % Refresh sensor bar charts (if requested):
            if this.showSensor
                set(this.hSensors,{'XData','YData'},{[mesh.elements.x],double([mesh.elements.isTroubled])})
            end
            % Refresh limiter bar charts (if requested):
            if this.showLimiter
                set(this.hLimiters,'XData',[mesh.elements.x])
                set(this.hLimiters,{'YData'},reshape(permute(num2cell(this.getLimiterData(mesh,this.rows),1),[3 2 1]),[],1))
            end
            % Refresh solution line plots:
            for k = 1:mesh.elementCount
                % Coordinates:
                x1 = mesh.elements(k).getNodeCoords';
                set(this.hNodes(:,k),'XData',x1)
                x2 = mesh.elements(k).getBreakCoords;
                set(this.hBreaks(:,k),'XData',x2)
                x3 = linspace(mesh.elements(k).xL,mesh.elements(k).xR,mesh.elements(k).dofCount*this.nPerDof);
                x = unique([x3 x2 x1],'sorted');
                set(this.hDiscrete(:,k),'XData',x)
                set(this.hExact(:,k),'XData',x)
                % Samples:
                j = cellfun(@(b) arrayfun(@(a) find(x == a,1),b),{x1 x2 x3},'UniformOutput',false);
                y = mesh.elements(k).interpolateStateAtCoords(x);
                set(this.hNodes(:,k),{'YData'},num2cell(y(this.rows,j{1}),2))
                set(this.hBreaks(:,k),{'YData'},num2cell(y(this.rows,j{2}),2))
                set(this.hDiscrete(:,k),{'YData'},num2cell(y(this.rows,:),2))
                y = this.solver.exactSolution(this.solver.timeNow,x);
                set(this.hExact(:,k),{'YData'},num2cell(y(this.rows,:),2))
            end
            % Append new norm samples to animated line plots:
            if ~isempty(this.norms)
                this.norms.compute(mesh,@(x) this.solver.exactSolution(this.solver.timeNow,x));
                for j = 1:numel(this.norms)
                    for i = this.rows
                        addpoints(this.hNorms(i,j),this.solver.timeNow,this.norms(j).vals(i));
                    end
                end
            end
            % Redraw:
            %%%this.updateYLims
            drawnow limitrate
        end
    end
    methods (Access = protected)
        %% Vertical axis limits
        function updateYLims(this)
            % Sets this Monitor's axis limits such that their range is 
            % allowed to increase but not shrink as the solver progresses.
            %
            % Retrieve updated solution values (exact and approximate):
            y = [...
                reshape([this.hExact.YData],numel(this.rows),[])...
                reshape([this.hDiscrete.YData],numel(this.rows),[])...
                ];
            % Non-shrinking axis limits update:
            this.ylims = [min([y this.ylims(:,1)],[],2) max([y this.ylims(:,2)],[],2)];
            set(this.hAxes(:,1),{'YLim'},num2cell(this.ylims,2))
        end
    end
    methods (Access = protected, Static)
        %% Process limiter data
        function kji = getLimiterData(mesh,eqs)
            % Processes th isLimited fields of a given mesh to facilitate
            % plotting limiter activation status as a stacked bar chart.
            % 
            % Preallocate and pad a 3D array of limiter activation flags:
            kji = zeros(mesh.elementCount,mesh.maxBasisCount-1,numel(eqs));
            % Retrieve N_k-1 highest mode flags (in descending order):
            for k = 1:mesh.elementCount
                kji(k,1:mesh.elements(k).dofCount-1,1:numel(eqs)) =...
                    permute(mesh.elements(k).isLimited(eqs,end:-1:2),[3 2 1]);
            end
            % Normalize by each element's maximum number of limited modes:
            kji = kji./([mesh.elements.dofCount]'-1);
        end
    end
end