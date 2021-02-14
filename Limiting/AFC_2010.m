classdef AFC_2010 < Limiter
    % Sub-cell a priori limiting for DG bases that satisfy the "partition 
    % of unity" requirement (e.g. Bspline). AFC approach to linearized FCT,
    % a la Kuzmin et al., 2010. Includes prelimiting and failsafe. Applied
    % to control variables. Full inter-patch coupling.
    %
    % For this version of AFC, the isLimited property is ill-defined, so it
    % is left set to false.
    %
    properties
        failsafeStages
        controlVars
    end
    properties (Access = protected)
        % Synchronizing functions:
        syncStatesFun = @AFC_2010.sync_skip
        syncFluxesFun = @AFC_2010.sync_skip
        invSyncStatesFun = @AFC_2010.sync_skip
        invSyncFluxesFun = @AFC_2010.sync_skip
    end
    methods
        %% Constructor
        function this = AFC_2010(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'FailsafeStages',0,@(x) validateattributes(x,{'double'},{'integer','scalar'}));
            addParameter(p,'ControlVars',[],@(x) validateattributes(x,{'double'},{'integer'}));
            % Parse the input sensor:
            parse(p,varargin{:});
            this.failsafeStages = p.Results.FailsafeStages;
            this.controlVars = p.Results.ControlVars;
        end
        %% Apply (each stage, override)
        function applyStage(~,varargin)
            % Intentionally does nothing.
        end
        %% Apply (full time-step, override)
        function applyStep(this,mesh,solver)
            % Corrects a low-order predictor (assumes AFC-supported basis
            % has been used in each mesh element).
            %
            mesh.computeResiduals(this.physics,solver) % extra residual evaluation (linearized high-order predictor)
            this.computeAntidiffusiveFluxes(mesh.elements,solver.timeDelta)
            this.applyAFC(mesh,solver) % linearised AFC
        end
        %% Apply (initialization, override)
        function applyInitial(this,mesh,solver)
            % Corrects a lumped (low-order) projection to make it 
            % constrained (high resolution). Same as step-wise (see above).
            %
            % Physics:
            this.physics = solver.physics;
            % Control variables:
            if isempty(this.controlVars)
                this.controlVars = 1:this.physics.equationCount;
            end
            % Sync functions:
            if isa(this.physics,'Euler')
                this.syncStatesFun = @AFC_2010.syncStates_euler;
                this.syncFluxesFun = @AFC_2010.syncFluxes_euler;
                this.invSyncStatesFun = @AFC_2010.invSyncStates_euler;
                this.invSyncFluxesFun = @AFC_2010.invSyncFluxes_euler;
            end
            % Constrained initialization:
            this.computeAntidiffusiveFluxes(mesh.elements,1)
            this.applyAFC(mesh,solver)
        end
        %% Name (extension)
        function name = getName(this)
            name = strrep(this.getName@Limiter,'_2010','');
            name = sprintf('%s, vars.: [%s]',name,num2str(this.controlVars));
            if this.failsafeStages > 0
                name = sprintf('%s + %d-stage failsafe',name,this.failsafeStages);
            end
        end
    end
    methods (Static, Access = {?AFC_2010,?EulerFCT})
        %% Compute antidiffusive fluxes
        function computeAntidiffusiveFluxes(elements,timeDelta)
            % Function that evaluates the antidiffusive fluxes of an array
            % of elements. They are stored as a sparse 2D array with each 
            % column being associated to a control point pair. All non-edge
            % pairs are padded with zeros.
            %
            % These fluxes are normalized by the jacobian of the mapping to
            % reference element space; this means that they are to be
            % divided by the lumped mass matrix IN REFERENCE SPACE.
            %
            for element = elements
                if ~isa(element.basis,'DGIGA_AFC')
                    error('AFC limiting is not supported for this basis.')
                end
                for edge = element.basis.edges
                    r = edge(1);
                    j = edge(2);
                    rj = edge(3);
                    element.antidiffusiveFluxes(:,rj) = timeDelta*(...
                        element.basis.massMatrix(r,j).*(element.residuals(:,r) - element.residuals(:,j)) +...
                        2/element.dx*element.diffusions{r,j}*(element.states(:,r) - element.states(:,j))...
                        );
                end
            end
        end
    end
    methods (Static, Access = protected)
        %% Zalesak's algorithm (vectorized)
        function alphas = zalesak(masses,fluxes,states,maxima,minima)
            % Generalized Zalesak's algorithm for the computation of a 
            % synchronized correction factor corresponding to a single 
            % control variable. From Kuzmin et al. 2010 (eqs. 8-11). Does
            % not include prelimiting.
            %
            % Instead of unrolling the loops over edges, applies fully
            % vectorized statements to a zero-padded 2D (sparse) matrix of
            % antidiffusive fluxes.
            %
            % Arguments
            %  masses: row array of lumped mass matrix diagonal entries (column: control point)
            %  fluxes: (sparse) row array of antidiffusive fluxes (column: edge; padded with zeros)
            %  states: row array of predictor state vectors (column: control point)
            %  maxima: row array of control point local maxima (column: control point)
            %  minima: idem, for minima
            %
            % Output
            %  alphas: (full) row array of FCT correction factors (column: edge)
            %
            % Set up:
            N = length(masses);
            fluxes = reshape(fluxes,N,N); % row: recieving control point (r); column: giving control point (j)
            % Net antidiffusive fluxes (1st column: +; 2nd column: -):
            P = [sum(max(0,fluxes),2) sum(min(0,fluxes),2)];
            % Distances to local extrema (1st column: +; 2nd column: -):
            Q = ([maxima; minima] - states)';
            % Control point correction factors (1st column: +; 2nd column: -):
            R = Q./P;
            R(~P) = inf; % all -Inf should now become +Inf
            R = min(1,masses'.*R);
            % Flux correction factors:
            alphas = R(:,1).*(fluxes >= 0) + R(:,2).*(fluxes < 0); %%% >= is safe iff prelimiting conservative fluxes (14/02/2020)
            alphas = min(alphas,alphas');
            % Re-arrange to match conservative antidiffusive flux arrangement:
            alphas = reshape(full(alphas),1,N^2);
        end
        %% Default sync function
        function sync_skip(~,~)
            % Placeholder method that intentionally does no conversion
            % to/from an element's staes and/or antidiffusive fluxes.
            % To be used for scalar physics.
            %
            % Intentionally does nothing.
        end
        function syncFluxes_prelimiting(element)
            % Does no conversion but applies prelimiting to the 
            % (conservative) antidiffusive fluxes present in the given 
            % element.
            %
            [I,N] = size(element.states);
            for i = I
                fluxes = reshape(element.antidiffusiveFluxes(i,:),N,N);
                fluxes(fluxes.*(element.states(i,:)-element.states(i,:)') > 0) = 0;
                element.antidiffusiveFluxes(i,:) = reshape(fluxes,1,N^2);
            end
        end
        %% Euler sync function (states)
        function syncStates_euler(element)
            % Converts the state vectors of a given element from
            % conservative variables (density, momentum, total energy) to
            % primitive variables (density, velocity, pressure), for the
            % case of 1D Euler's equations.
            %
            element.states = Euler.stateToPrimitive(element.states);
        end
        %% Euler inverse sync function (states)
        function invSyncStates_euler(element)
            % Converts the state vectors of each element from primitive to
            % conservative variables, for the case of 1D Euler's equations.
            %
            element.states = Euler.primitiveToState(element.states);
        end
        %% Euler sync function (fluxes)
        function syncFluxes_euler(element)
            % Converts the antidiffusive flux vectors of a given element 
            % from conservative variable differences (density, momentum, 
            % total energy) to primitive ones (density, velocity,
            % pressure), for the case of 1D Euler's equations. A la Kuzmin
            % et al. 2010, eq. 6.
            %
            % Assumes that state vectors contain primitive variables.
            %
            r = element.basis.edges(1,:);
            rj = element.basis.edges(3,:);
            % Pressure flux:
            element.antidiffusiveFluxes(3,rj) = .4*(...
                element.antidiffusiveFluxes(3,rj) +...
                .5*element.states(2,r).^2.*element.antidiffusiveFluxes(1,rj) -...
                element.states(2,r).*element.antidiffusiveFluxes(2,rj));
            % Velocity flux:
            element.antidiffusiveFluxes(2,rj) = (...
                element.antidiffusiveFluxes(2,rj) -...
                element.states(2,r).*element.antidiffusiveFluxes(1,rj)...
                )./element.states(1,r);
        end
        %% Euler inverse sync function (fluxes)
        function invSyncFluxes_euler(element)
            % Converts the antidiffusive flux vectors of a given element 
            % from primitive variable differences to conservative ones, for
            % the 1D Euler's equations.
            %
            % Assumes that state vectors contain primitive variables.
            %
            r = element.basis.edges(1,:);
            rj = element.basis.edges(3,:);
            % Momentum flux:
            element.antidiffusiveFluxes(2,rj) =...
                element.antidiffusiveFluxes(2,rj).*element.states(1,r) +...
                element.states(2,r).*element.antidiffusiveFluxes(1,rj);
            % Total energy flux:
            element.antidiffusiveFluxes(3,rj) =...
                2.5*element.antidiffusiveFluxes(3,rj) -...
                .5*element.states(2,r).^2.*element.antidiffusiveFluxes(1,rj) +...
                element.states(2,r).*element.antidiffusiveFluxes(2,rj);
        end
    end
    methods (Access = protected)
        %% Apply (AFC)
        function applyAFC(this,mesh,solver)
            % Employs AFC a priori limiting a la Kuzmin et al., 2010 on a
            % "quasi-nodal" DGIGA discretization. Interpatch coupling is 
            % made via "stitching" of basis function supports across 
            % breakpoints.
            %
            % Assumes that the mesh contains the low-order predictor at the
            % next time-step ("transported and diffused") as well as some 
            % antidiffusive fluxes to constrain. It then reconstructs a 
            % limited high order solution (which is LED) via an AFC-based,
            % sequential, synchronized, linearized FCT procedure applied
            % to the antidiffusive flux components associated to each
            % combination of two control points. Additionally, a failsafe
            % check is carried out afterwards.
            %
            % Apply sensor to (linearized) high-order solutions:
            this.applySensor(mesh,solver)
            % Retrieve troubled elements:
            isTroubled = [mesh.elements.isTroubled];
            elements = mesh.elements(isTroubled(:,:,this.priority));
            % Prelimiting (in conserved variables):
            this.applyPrelimiting(elements)
            % Determine local extrema (in control variables):
            this.findExtrema(mesh)
            % FCT limiting (in control variables):
            this.applySynchronizedFCT(elements)
            % Failsafe limiting (in control variables):
            this.applyFailsafe(elements)
        end
        %% Apply sensor
        function applySensor(this,mesh,solver)
            % The predictor-corrector nature of AFC makes it necessary to 
            % apply the sensor as follows. It is assumed that antidiffusive
            % fluxes and low-order predictors exist in each element of the
            % given mesh.
            %
            % Apply raw antidiffusive fluxes to each element:
            for element = mesh.elements
                element.applyAntidiffusiveFluxes;
                element.isLimited(:,:,this.priority) = false(size(element.states));
            end
            % Detect troubled cells in the linearised high-order solution:
            this.sensor.apply(mesh,solver,this.priority)
            % Recover low-order predictors of troubled patches:
            isTroubled = [mesh.elements.isTroubled];
            for element = mesh.elements(isTroubled(:,:,this.priority))
                element.removeAntidiffusiveFluxes(1);
            end
        end
        %% Prelimiting, conservative variables
        function applyPrelimiting(this,elements)
            % Prelimiting (a la Kuzmin 2012, eq. 79). Applied to
            % antidiffusive fluxes of conservative variables.
            %
            for element = elements
                N = element.basis.basisCount;
                M = N^2;
                % Loop over PDE components:
                for i = this.physics.equationCount
                    fluxes = reshape(element.antidiffusiveFluxes(i,:),N,N);
                    fluxes(fluxes.*(element.states(i,:) - element.states(i,:)') > 0) = 0;
                    element.antidiffusiveFluxes(i,:) = reshape(fluxes,1,M);
                end
            end
        end
        %% Deduce local extrema (full support across patch boundaries)
        function findExtrema(this,mesh)
            % This method initializes arrays maxima/minima of each troubled
            % element with its local maximum/minimum control variable
            % values. Only troubled elements have their conserved variables
            % replaced with control variables; non-troubled neighbors are
            % returned to conserved variables.
            %
            % Searches extrema in the full stencil width across each patch 
            % interface, i.e. all basis modes that "touch" a patch 
            % interface (on both sides of it) share common extrema.
            %
            % Find intra-patch extrema (also ghost elements):
            for element = [mesh.elements mesh.boundaries.ghostElement]
                % Transform to control variables:
                this.syncStatesFun(element)
                % Aliases:
                r = element.basis.pairs(1,:); % test functions
                j = element.basis.pairs(2,:); % overlapping basis functions
                % Loop over control points of current patch:
                for m = r
                    ids = j(r == m);
                    element.maxima(:,m) = max(element.states(:,ids),[],2);
                    element.minima(:,m) = min(element.states(:,ids),[],2);
                end
            end
            % Communicate inter-patch extrema to troubled elements:
            %
            %  Right(left)-most control point extrema of patch k is also
            %  the extrema of the first(last) control point of patch k+1
            %  (k-1).
            %
            for edge = mesh.edges % loop over edges
                % Common maxima:
                aux = max([edge.elementL.maxima(:,end),edge.elementR.maxima(:,1)],[],2);
                edge.elementL.maxima(:,end) = aux;
                edge.elementR.maxima(:,1) = aux;
                % Common minima:
                aux = min([edge.elementL.minima(:,end),edge.elementR.minima(:,1)],[],2);
                edge.elementL.minima(:,end) = aux;
                edge.elementR.minima(:,1) = aux;
            end
            % Distribute inter-patch extrema inwards of each troubled cell:
            isTroubled = [mesh.elements.isTroubled];
            mask = isTroubled(:,:,this.priority);
            for element = mesh.elements(mask)
                % Aliases:
                basis = element.basis;
                r = basis.edges(1,:); % test functions
                j = basis.edges(2,:); % overlapping basis functions (excluding self-support)
                % Left/right-most control point nearest neighbours:
                for m = [1 basis.basisCount]
                    idsTo = j(r == m);
                    idsFrom = m*ones(size(idsTo));
                    element.maxima(:,idsTo) = max(element.maxima(:,idsTo),element.maxima(:,idsFrom));
                    element.minima(:,idsTo) = min(element.minima(:,idsTo),element.minima(:,idsFrom));
                end
            end
            % Return untroubled cells to conservative variables:
            for element = mesh.elements(~mask)
                this.invSyncStatesFun(element)
            end
            % Override extrema of control points closest to mesh boundaries
            % (Kuzmin et al, 2012; remark 5, pp. 163, bottom):
%             mesh.elements(1).maxima(:,1) = inf;
%             mesh.elements(1).minima(:,1) = -inf;
%             mesh.elements(end).maxima(:,end) = inf;
%             mesh.elements(end).minima(:,end) = -inf;
        end
        %% Linearized FCT, synchronized
        function applySynchronizedFCT(this,elements)
            % Linearized FCT algorithm based on Kuzmin et al. 2010. 
            % Synchronized correction using all primitive variables as
            % control variables.
            %
            % Assumes that each element in the mesh contains state and 
            % residual vectors (2D arrays), as well as AFC diffusion matrix
            % blocks (cell array) and antidiffusive fluxes (2D array padded 
            % with zeros), all of them associated with the low-order 
            % predictor at the next time-level (transported and diffused). 
            % 
            % After execution, the mesh contains the corrected state 
            % vectors (in conservative variables).
            %
            % Control variable loop:
            for i = this.controlVars
                % Patch loop:
                for element = elements
                    % Convert antidiffusive fluxes to control variables:
                    this.syncFluxesFun(element)
                    % Zalesak's scalar algorithm (single control variable):
                    alphas = AFC_2010.zalesak(...
                        element.basis.lumpedMassMatrixDiagonal,...
                        element.antidiffusiveFluxes(i,:),...
                        element.states(i,:),...
                        element.maxima(i,:),...
                        element.minima(i,:));
                    % Revert to conservative antidiffusive fluxes:
                    this.invSyncFluxesFun(element)
                    % Limit conservative antidiffusive fluxes:
                    element.antidiffusiveFluxes = alphas.*element.antidiffusiveFluxes;
                    %%%element.isLimited(i,:,this.priority) = sum(reshape(1-alphas,element.dofCount,element.dofCount),2)';
                end
            end
            % Finalize:
            for element = elements
                % Revert state to conservative variables:
                this.invSyncStatesFun(element)
                % Apply limited antidiffusive fluxes to each element:
                element.applyAntidiffusiveFluxes
            end
        end
        %% Failsafe flux limiting
        function applyFailsafe(this,elements)
            % Failsafe flux correction based on local imposition of FCT
            % constraints to control variables, a la Kuzmin et al. 2010.
            %
            % Assumes that each element in the mesh contains state vectors
            % (in conservative variables; 2D array) associated with a
            % low-order predictor at the next time-level. After execution,
            % the mesh contains the corrected solution states
            % (conservative variables).
            %
            % Limiting factors 'beta' are computed sequentially, e.g.:
            % first, the density one is obtained and used to remove
            % part of the raw (conservative) antidiffusive fluxes; the
            % resulting "partially corrected solution" is then used to
            % calculate a pressure limiting factor, and so on.
            %
            % Loop over stages:
            for m = 1:this.failsafeStages
                % Loop over elements:
                for element = elements
                    % Preallocate the beta factors into a 2D array, such that
                    % each row represents control point 'r', and each column
                    % control point 'j'.
                    N = element.basis.basisCount;
                    betas = zeros(N,N);
                    % Loop over control variables:
                    for i = this.controlVars
                        % Convert state to control variables:
                        this.syncStatesFun(element);
                        % Compute failsafe limiting coefficents associated
                        % with the current control variable:
                        isInvalid = element.states(i,:) < element.minima(i,:) - 1e-6 | element.states(i,:) > element.maxima(i,:) + 1e-6;
                        betas(isInvalid | isInvalid') = m/this.failsafeStages;
                        % Revert states back to conservative variables:
                        this.invSyncStatesFun(element);
                        % Substract a fraction of antidiffusive fluxes:
                        element.removeAntidiffusiveFluxes(betas);
                    end
                end
            end
        end
    end
end