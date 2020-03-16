classdef AFC < Limiter
    properties
        everyStage = false; % false -> too diffusive; true -> too antidiffusive
        physics
        % Synchronizing functions:
        syncStatesFun = @AFC.sync_skip
        syncFluxesFun = @AFC.sync_skip %%% AFC.syncFluxes_prelimiting
        invSyncStatesFun = @AFC.sync_skip
        invSyncFluxesFun = @AFC.sync_skip
    end
    methods
        %% Constructor
        function this = AFC(varargin)
            this.physics = physics;
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'Sensor',Sensor,@(x) validateattributes(x,{'Sensor'},{}));
            % Parse the input sensor:
            parse(p,varargin{:});
            this.sensor = p.Results.Sensor;
        end
        %% Algebraic Flux Correction (linearized FCT + failsafe)
        function apply(this,mesh,solver,isInitial)
            % Employs AFC-based FCT a priori limiting a la Kuzmin et al.,
            % 2010 on a "quasi-nodal" DGIGA discretization. Interpatch
            % coupling is made via P and Q communication across edges.
            %
            % Assumes that the mesh contains the low order solution at the
            % next time-step ("transported and diffused"). It then 
            % reconstructs a limited high order solution (which is LED) via
            % an AFC-based, sequential, synchronized FCT procedure applied
            % to the antidiffusive flux components associated to each
            % combination of two control points. Additionally, a failsafe
            % check is carried out.
            %
            % Valid for a system of conservation laws of any number of
            % components.
            %
            % Default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
            
            % Apply on every remaining element:
            for element = elements
                this.applyOnElement(element);
            end
            
            % Extra residual evaluation (future predictor estimate):
            mesh.computeResiduals(this.physics);
            % Compute raw antidiffusive fluxes (conservative variables):
            this.computeAntidiffusiveFluxes(mesh,timeDelta);
            % Apply prelimiting:
            this.applyPrelimiting(mesh);
            % Convert state to control variables:
            this.syncStatesFun(mesh);
            % Determine local predictor extrema:
            this.findExtrema(mesh);
            % Apply FCT limiting:
            this.applySynchronizedFCT(mesh);
            % Apply failsafe limiting:
            this.applyFailsafe(mesh);
        end
        %% Prelimiting, conservative variables
        function applyPrelimiting(this,mesh)
            % Prelimiting (a la Kuzmin 2012, eq. 79). Applied to
            % antidiffusive fluxes of conservative variables.
            %
            % Loop over elements:
            for element = mesh.elements
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
        %% Linearized FCT, synchronized
        function applySynchronizedFCT(this,mesh)
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
            for i = this.physics.controlVars
                % Patch loop:
                for element = mesh.elements
                    % Convert antidiffusive fluxes to control variables:
                    this.syncFluxesFun(element);
                    % Zalesak's scalar algorithm (single control variable):
                    alphas = AFC.zalesak(...
                        element.basis.lumpedMassMatrixDiagonal,...
                        element.antidiffusiveFluxes(i,:),...
                        element.states(i,:),...
                        element.maxima(i,:),...
                        element.minima(i,:));
                    % Revert to conservative antidiffusive fluxes:
                    this.invSyncFluxesFun(element);
                    % Limit conservative antidiffusive fluxes:
                    element.antidiffusiveFluxes = alphas.*element.antidiffusiveFluxes;
                end
            end
            % Revert state to conservative variables:
            this.invSyncStatesFun(mesh);
            % Apply limited antidiffusive fluxes to each element:
            for element = mesh.elements
                element.applyAntidiffusiveFluxes;
            end
        end
        %% Failsafe flux limiting
        function applyFailsafe(this,mesh,M)
            % Failsafe flux correction based on local imposition of FCT
            % constraints to control variables, a la Kuzmin et al. 2010.
            %
            % Assumes that each element in the mesh contains state vectors 
            % (in conservative variables; 2D array) associated with a 
            % low-order predictor at the next time-level. After execution,
            % the mesh contains the corrected solution states 
            % (conservative variables).
            %
            if nargin < 3
                M = 1; % default: single stage
            end
            % Loop over elements:
            for element = mesh.elements
                % Precompute some stuff:
                N = element.basis.basisCount;
                betas = sparse(N,N);
                % Apply the failsafe limiting in M stages:
                for m = 1:M
                    % Convert state to control variables:
                    this.syncStatesFun(mesh);
                    % Loop over control variables:
                    for i = 1:this.physics.equationCount
                        updateBetas(...
                            element.states(i,:),...
                            element.maxima(i,:),...
                            element.minima(i,:));
                    end
                    % Revert states back to conservative variables:
                    this.invSyncStatesFun(mesh);
                    % Substract a fraction of antidiffusive fluxes:
                    element.removeAntidiffusiveFluxes(betas);
                end
            end
            function updateBetas(states,maxima,minima)
                % Computes the fraction of antidiffusive fluxes to remove
                % in the current failsafe stage. Nested. Vectorized.
                %
                % Arguments:
                %  states: row array of predictor control variable states (column: control point)
                %  maxima: row array of control variable local maxima (column: control point)
                %  minima: idem, for minima
                %
                % Find indices of troubled control points:
                troubled = states < minima - 1e-6 | states > maxima + 1e-6;
                % Find indices of troubled antidiffusive fluxes:
                [r,j] = find(troubled.*troubled');
                % Update corrections of troubled antidiffusive fluxes:
                betas(r,j) = m/M;
            end
        end
    end
    methods (Static)
        %% Compute antidiffusive fluxes
        function computeAntidiffusiveFluxes(mesh,timeDelta)
            % Function that evaluates the antidiffusive fluxes from a given
            % mesh. They are stored element-wise as a sparse 2D array with
            % each column being associated to a control point pair. All
            % non-edge pairs are padded with zeros.
            %
            % These fluxes are normilized by the jacobian of the mapping to
            % reference element space; this means that they are to be
            % divided by the lumped mass matrix IN REFERENCE SPACE.
            %
            for element = mesh.elements
                basis = element.basis;
                for edge = basis.edges
                    r = edge(1);
                    j = edge(2);
                    rj = edge(3);
                    element.antidiffusiveFluxes(:,rj) = timeDelta*(...
                        basis.massMatrix(r,j).*(element.residuals(:,r) - element.residuals(:,j)) +...
                        2/element.dx*element.diffusions{r,j}*(element.states(:,r) - element.states(:,j))...
                        );
                end
            end
        end
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
        %% Deduce local extrema (full support across patch boundaries)
        function findExtrema(mesh)
            % This method initializes arrays maxima/minima of each element 
            % with its local maximum/minimum state values.
            %
            % Searches extrema in the full stencil width across each patch 
            % interface, i.e. all basis modes that "touch" a patch 
            % interface (on both sides of it) share common extrema.
            %
            % Find intra-patch extrema:
            N = mesh.elementCount;
            for element = mesh.elements
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
            % Communicate inter-patch extrema:
            %
            %  Right(left)-most control point extrema of patch k is also
            %  the extrema of the first(last) control point of patch k+1
            %  (k-1).
            %
            for k = 2:N % loop over (left) interfaces
                % Aliases:
                elementL = mesh.elements(k-1);
                elementR = mesh.elements(k);
                % Common maxima:
                aux = max([elementL.maxima(:,end),elementR.maxima(:,1)],[],2);
                elementL.maxima(:,end) = aux;
                elementR.maxima(:,1) = aux;
                % Common minima:
                aux = min([elementL.minima(:,end),elementR.minima(:,1)],[],2);
                elementL.minima(:,end) = aux;
                elementR.minima(:,1) = aux;
            end
            % Distribute inter-patch extrema inwards of each patch:
            for element = mesh.elements
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
            % Override extrema of control points closest to mesh boundaries:
            % (Kuzmin et al, 2012; remark 5, pp. 163, bottom)
            mesh.elements(1).maxima(:,1) = inf;
            mesh.elements(1).minima(:,1) = -inf;
            mesh.elements(end).maxima(:,end) = inf;
            mesh.elements(end).minima(:,end) = -inf;
        end
        %% Deduce local extrema (minimal support across patch boundaries)
        function findExtremaBAD(mesh)
            % This method initializes arrays maxima/minima of each element 
            % with its local maximum/minimum state values.
            %
            % Searches extrema only on the first control point across each 
            % patch interface, i.e. only the extrema of the right-most
            % mode of the left patch can be candidate extrema of the right
            % patch.
            %
            % Alias:
            N = mesh.elementCount;
            % Determine common extrema (nearest modes to patch boundaries):
            for k = 2:N % loop over (left) interfaces
                % Aliases:
                elementL = mesh.elements(k-1);
                elementR = mesh.elements(k);
                % Common maxima:
                aux = max([elementL.states(:,end),elementR.states(:,1)],[],2);
                elementL.maxima(:,N) = aux;
                elementR.maxima(:,1) = aux;
                % Common minima:
                aux = min([elementL.states(:,end),elementR.states(:,1)],[],2);
                elementL.minima(:,N) = aux;
                elementR.minima(:,1) = aux;
            end
            % Find intra-patch extrema:
            for element = mesh.elements
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
            % Distribute inter-patch extrema inwards of each patch:
            for element = mesh.elements
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
            % Override extrema of control points closest to mesh boundaries:
            % (Kuzmin et al, 2012; remark 5, pp. 163, bottom)
            mesh.elements(1).maxima(:,1) = inf;
            mesh.elements(1).minima(:,1) = -inf;
            mesh.elements(end).maxima(:,end) = inf;
            mesh.elements(end).minima(:,end) = -inf;
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
        function syncStates_euler(mesh)
            % Converts the state vectors of a given mesh from
            % conservative variables (density, momentum, total energy) to
            % primitive variables (density, velocity, pressure), for the
            % case of 1D Euler's equations.
            %
            for element = mesh.elements
                element.states = Euler.stateToPrimitive(element.states);
            end
        end
        %% Euler inverse sync function (states)
        function invSyncStates_euler(mesh)
            % Converts the state vectors of each element from primitive to
            % conservative variables, for the case of 1D Euler's equations.
            %
            for element = mesh.elements
                element.states = Euler.primitiveToState(element.states);
            end
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
end