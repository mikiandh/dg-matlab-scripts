classdef Monitor < handle
    % Class that handles monitoring tasks related to a particular solver
    % instance. These include plotting of the solution and computation of
    % solution and/or error norms.
    properties (Constant)
        colorSkip = {'w','k'} % colors to avoid in the colormaps
        nPerDof = 6 % number of sample points per DoF
        figureSizes = [700 400; 1400 800] % width, height (minimum); idem (maximum)
    end
    properties (SetAccess = immutable)
        solver
        norms
        rows
        cols
    end
    properties (Access = protected)
        % Colormaps:
        cDiscrete
        cLimiters
        cNorms
        % Handles to graphics objects:
        hFigure % handle to the figure
        hStaticTitle % handle to static title annotation
        hDynamicTitle % idem, dynamic (updated at each iteration)
        hPanel % handle to the UI panel
        hAxes % 2D array of axes (i.e. subplot) handles
        hDiscrete % 2D array of patch-wide approximate solutions (row: equation; column: cell)
        hExact % idem, exact solutions
        hControlPoints % idem, control point coordinates and values (if any)
        hNodes % idem, nodal coordinates and values
        hBreaks % idem, breakpoints
        hSensors % column array of sensor activation bar charts
        hLimiters % idem, limiter activation bar chars (stacked)
        hNorms % 2D array of instantaneous norm samples (row: equation; column: norm type)
        hLegend % handle to the legend of the norm subplots
        % Other graphics properties:
        figurePosition
        titleHeight = 80
        ylims = [0 1] % columns: min, max (#rows is set automatically)
    end
    methods
        %% Constructor
        function this = Monitor(solver,varargin)
            % Initialize parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'solver',@(x)validateattributes(x,{'Solver'},{}))
            addParameter(p,'norms',[],@(x)validateattributes(x,{'Norm'},{'vector'}))
            addParameter(p,'equations',1:solver.physics.equationCount,@(x)validateattributes(x,{'numeric'},{'vector','integer','>',0,'<=',solver.physics.equationCount}))
            p.parse(solver,varargin{:});
            % Set properties from parser:
            this.solver = p.Results.solver;
            this.norms = unique(p.Results.norms);
            % Subplot indices:
            this.rows = p.Results.equations;
            this.cols = 1:1 + ~isempty(this.norms);
            % Repmat initial axis limits:
            this.ylims = repmat(this.ylims,numel(this.rows),1);
            % Set main figure position and size:
            this.figurePosition = [400 100 min(this.figureSizes.*[this.cols(end) this.rows(end);1 1])];
            % Normalize title height:
            this.titleHeight = this.titleHeight/this.figurePosition(4);
        end
        %% Initialize
        function initialize(this,mesh,fun)
            % Initializes a figure with a 2D arrangement of subplots.
            %
            % The left column shows one solution component per row, over
            % the entire domain, at a particular instant in time. Each
            % patch (i.e. DG element, a.k.a macro-element) has a different
            % colour. Breakpoints are indicated by vertical cross markers.
            % Nodes (if any) are shown as dot markers. Control points (if
            % any) by diagonal cross markers.
            %
            % The right column shows samples over time of component-wise 
            % solution and/or error norms, over the entire domain. Each
            % norm has a different colour. Sample points are indicated by
            % cross markers.
            %
            % Initialize a figure:
            this.hFigure = figure('Position',this.figurePosition,'Name',class(this));
            % Instantiate a uipanel:
            this.hPanel = uipanel(this.hFigure,'Position',[0 0 1 1-this.titleHeight],'BackgroundColor','w','BorderType','none');
            % Initialize "supertitle" annotations:
            this.hStaticTitle = annotation(this.hFigure,'textbox',[0 1-this.titleHeight 1 this.titleHeight],...
                'String',{[this.solver.physics.getInfo '; ' mesh.boundaries.getInfo],[mesh.getInfo '; ' this.solver.limiter.getInfo]},...
                'FontWeight','bold','EdgeColor','none',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            this.hDynamicTitle = annotation(this.hFigure,'textbox',[0 1-this.titleHeight 1 this.titleHeight],...
                'String',this.solver.getInfo,...
                'FontWeight','bold','EdgeColor','none',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
            % Preallocate graphics objects:
            this.hDiscrete = gobjects(this.rows(end),mesh.elementCount);
            this.hExact = this.hDiscrete;
            this.hControlPoints = this.hDiscrete;
            this.hNodes = this.hDiscrete;
            this.hBreaks = this.hDiscrete;
            this.hSensors = gobjects(this.rows(end),1);
            this.hLimiters = gobjects(this.rows(end),mesh.maxBasisCount);
            this.hNorms = gobjects(this.rows(end),numel(this.norms));
            % Initialize 2D lattice of subplot axes:
            k = 0;
            for i = this.rows
                for j = this.cols
                    k = k + 1;
                    this.hAxes(i,j) = subplot(this.rows(end),this.cols(end),k,'Parent',this.hPanel);
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
            this.cLimiters = distinguishable_colors(mesh.maxBasisCount,this.colorSkip);
            this.cDiscrete = distinguishable_colors(mesh.elementCount,this.colorSkip);
            this.cNorms = distinguishable_colors(numel(this.norms),this.colorSkip);
            % Initialize sensor bar charts:
            for i = this.rows
                this.hSensors(i) = bar(this.hAxes(i,1),[mesh.elements.x],[mesh.elements.isTroubled],1);
            end
            set(this.hSensors,'FaceAlpha',.05,'EdgeAlpha',.1,...
                'EdgeColor',this.cLimiters(1,:),'FaceColor',this.cLimiters(1,:))
            if strcmp(this.solver.limiter.sensor,'Sensor')
                [this.hSensors.Visible] = deal('off');
            end
            % Initialize limiter bar charts:
            for i = this.rows
                this.hLimiters(i,:) = bar(this.hAxes(i,1),this.getLimiterData(mesh,i),1,'stacked'); % fraction of 1
                set(this.hLimiters(i,:),'XData',[mesh.elements.x]);  % hotfix around the "equal lengths bug"
                set(this.hLimiters(i,:),{'EdgeColor'},num2cell(this.cLimiters,2),{'FaceColor'},num2cell(this.cLimiters,2));
            end
            set(this.hLimiters,'FaceAlpha',.1,'EdgeAlpha',.2)
            if strcmp(this.solver.limiter,'Limiter')
                [this.hSensors.Visible] = deal('off');
            end
            % Initialize solution vs. space plots:
            for k = 1:mesh.elementCount
                % Sample coordinates:
                x0 = mesh.elements(k).getControlCoords;
                x1 = mesh.elements(k).getNodeCoords';
                x2 = mesh.elements(k).getBreakCoords;
                x3 = linspace(mesh.elements(k).xL,mesh.elements(k).xR,ceil(mean(1000,mesh.elements(k).dofCount*this.nPerDof)));
                x = unique([x3 x2 x1],'sorted');
                % Indices of sample subsets:
                j = cellfun(@(b) arrayfun(@(a) find(x == a,1),b),{x1 x2 x3},'UniformOutput',false);
                % Solution samples:
                t = mesh.elements(k).states(:,1:numel(x0));
                y = mesh.elements(k).interpolateStateAtCoords(x);
                z = fun(x); % passed by the solver (either initial condition or exact solution)
                for i = this.rows
                    this.hExact(i,k) = plot(this.hAxes(i,1),x,z(i,:),'LineStyle','--','Color','k');
                    this.hDiscrete(i,k) = plot(this.hAxes(i,1),x,y(i,:),'LineStyle','-','Color',this.cDiscrete(k,:));
                    this.hBreaks(i,k) = plot(this.hAxes(i,1),x2,y(i,j{2}),'LineStyle','none','Marker','+','Color',this.cDiscrete(k,:));
                    this.hNodes(i,k) = plot(this.hAxes(i,1),[nan x1],[nan y(i,j{1})],'LineStyle','none','Marker','.','Color',this.cDiscrete(k,:));
                    this.hControlPoints(i,k) = plot(this.hAxes(i,1),[nan x0],[nan t(i,:)],'LineStyle','none','Marker','x','Color',this.cDiscrete(k,:));
                end
            end
            % Hide control points by default:
            [this.hControlPoints.Visible] = deal('off');
            % Initialize norm vs. time plots (if any):
            if ~isempty(this.norms)
                this.norms.compute(mesh,@(x) this.solver.exactSolution(this.solver.timeNow,x));
                for j = 1:numel(this.norms)
                    for i = this.rows
                        this.hNorms(i,j) = animatedline(this.solver.timeNow,this.norms(j).vals(i),'Parent',this.hAxes(i,2),...
                            'LineStyle','-','Marker','x','Color',this.cNorms(j,:),'DisplayName',char(this.norms(j)));
                    end
                end
                % Initialize legend:
                this.hLegend = legend(this.hAxes(end),'Location','NorthEast');
                this.hLegend.ItemHitFcn = @this.toggle_norms;
            end
            % Set axis limits:
            set(this.hAxes(:,1),'XLim',[mesh.edges(1).coord mesh.edges(end).coord])
            try
                set(this.hAxes(:,2),'XLim',[this.solver.timeSpan this.solver.timeStop])
            catch
                % Nothing serious. Keep going.
            end
            %%%this.updateYLims
            % Initialize a customized toolbar:
            this.setToolbar
            % Redraw:
            drawnow
        end
        %% Refresh
        function refresh(this,mesh)
            % Replots the monitor figure with current solution and norms.
            %
            % Refresh the (dynamic) title:
            this.hDynamicTitle.String = this.solver.getInfo;
            % Refresh sensor bar charts:
            set(this.hSensors,{'XData','YData'},{[mesh.elements.x],double([mesh.elements.isTroubled])})
            % Refresh limiter bar charts:
            set(this.hLimiters,'XData',[mesh.elements.x])
            set(this.hLimiters,{'YData'},reshape(permute(num2cell(this.getLimiterData(mesh,this.rows),1),[3 2 1]),[],1))
            % Refresh solution line plots:
            for k = 1:mesh.elementCount
                % Coordinates:
                x0 = mesh.elements(k).getControlCoords;
                set(this.hControlPoints(:,k),'XData',x0)
                x1 = mesh.elements(k).getNodeCoords';
                set(this.hNodes(:,k),'XData',x1)
                x2 = mesh.elements(k).getBreakCoords;
                set(this.hBreaks(:,k),'XData',x2)
                x3 = linspace(mesh.elements(k).xL,mesh.elements(k).xR,ceil(mean(1000,mesh.elements(k).dofCount*this.nPerDof)));
                x = unique([x3 x2 x1],'sorted');
                set(this.hDiscrete(:,k),'XData',x)
                set(this.hExact(:,k),'XData',x)
                % Samples:
                set(this.hControlPoints(:,k),{'YData'},num2cell(mesh.elements(k).states(this.rows,1:numel(x0)),2));
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
        %% Initialize custom toolbar
        function setToolbar(this,varargin)
            % Add this Monitor's custom tools into its figure's toolbar.
            % Retrive existing toolbar:
            hToolbar = findall(this.hFigure,'type','uitoolbar');
            % Load icons:
            icons = load('monitor_icons.mat','-regexp','icon');
            % Toggle visibility tools:            
            uitoggletool(hToolbar,'State',this.hDiscrete(1).Visible,'ClickedCallback',@this.toggle_discrete,'Tooltip','Show discrete solution','CData',icons.icon_discrete,'Separator','on')
            uitoggletool(hToolbar,'State',this.hExact(1).Visible,'ClickedCallback',@this.toggle_exact,'Tooltip','Show exact solution','CData',icons.icon_exact)
            uitoggletool(hToolbar,'State',this.hControlPoints(1).Visible,'ClickedCallback',@this.toggle_controlPoints,'Tooltip','Show control points','CData',icons.icon_controlPoints)
            uitoggletool(hToolbar,'State',this.hNodes(1).Visible,'ClickedCallback',@this.toggle_nodes,'Tooltip','Show nodes','CData',icons.icon_nodes)
            uitoggletool(hToolbar,'State',this.hBreaks(1).Visible,'ClickedCallback',@this.toggle_breaks,'Tooltip','Show breakpoints','CData',icons.icon_breaks)
            uitoggletool(hToolbar,'State',this.hSensors(1).Visible,'ClickedCallback',@this.toggle_sensor,'Tooltip','Show sensor','CData',icons.icon_sensor)
            uitoggletool(hToolbar,'State',this.hLimiters(1).Visible,'ClickedCallback',@this.toggle_limiter,'Tooltip','Show limiter','CData',icons.icon_limiter)
        end
        %% Toggle norm visibility
        function toggle_norms(this,~,event)
            % Callback function that toggles visibility of all animated
            % lines corresponding to the norm that has just been clicked on
            % by the user, on the legend.
            [~,j] = find(this.hNorms == event.Peer);
            if strcmp(event.Peer.Visible,'on')
                [this.hNorms(:,j).Visible] = deal('off');
            else
                [this.hNorms(:,j).Visible] = deal('on');
            end
        end
        %% Toggle breakpoint visibility
        function toggle_breaks(this,src,~)
            % Callback function that switches breakpoint visibility.
            [this.hBreaks.Visible] = deal(src.State);
        end
        %% Toggle node visibility
        function toggle_nodes(this,src,~)
            % Callback function that switches node visibility.
            [this.hNodes.Visible] = deal(src.State);
        end
        %% Toggle control point visibility
        function toggle_controlPoints(this,src,~)
            % Callback function that switches control point visibility.
            [this.hControlPoints.Visible] = deal(src.State);
        end
        %% Toggle discrete solution visibility
        function toggle_discrete(this,src,~)
            % Callback function that switches visibility of the piecewise
            % polynomial approximation to the solution.
            [this.hDiscrete.Visible] = deal(src.State);
        end
        %% Toggle exact solution visibility
        function toggle_exact(this,src,~)
            % Callback function that switches exact solution visibility.
            [this.hExact.Visible] = deal(src.State);
        end
        %% Toggle sensor visibility
        function toggle_sensor(this,src,~)
            % Callback function that switches sensor visibility.
            [this.hSensors.Visible] = deal(src.State);
        end
        %% Toggle limiter visibility
        function toggle_limiter(this,src,~)
            % Callback function that switches limiter visibility.
            [this.hLimiters.Visible] = deal(src.State);
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
                kji(k,1:mesh.elements(k).dofCount,1:numel(eqs)) =...
                    permute(mesh.elements(k).isLimited(eqs,end:-1:1),[3 2 1]);
            end
            % Normalize by each element's maximum number of limited modes:
            kji = kji./([mesh.elements.dofCount]'-1);
        end
    end
end