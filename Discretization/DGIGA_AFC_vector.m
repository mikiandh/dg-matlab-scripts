classdef DGIGA_AFC_vector < Bspline
    properties
        lumpedMassMatrixDiagonal % non-zero entries of the lumped mass matrix
        modePairs % 2D array of test function/basis function pairs within mutual nonzero support; 1st row: test function, r; 2nd row: basis function, j
        edges % idem, but excludes pairs of repeated control points ("virtual" edges a la Kuzmin et al. 2012)
    end
    methods
        %% Constructor
        function this = DGIGA_AFC_vector(varargin)
            this@Bspline(varargin{:});
            if nargin > 0
                this.lumpedMassMatrixDiagonal = full(sum(this.massMatrix,1));
                [r,j] = find(this.massMatrix);
                this.modePairs = vertcat(r',j');
                this.edges = this.modePairs(:,r~=j);
            end
        end
        %% Instantiate from prototype
        function this = clone(prototype,degree)
            this = DGIGA_AFC_vector(prototype.knots,degree,prototype.smoothness);
        end
        %% DGIGA-AFC (low order predictor) operator
        function computeResiduals(this,element,physics)
            
%%% NOT WORKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(element.jacobians)
                  this.initializeBlockMatrices(element,this.basisCount,physics.equationCount); % TO DO: place somewhere else (projection?)
            end
            this.computeDiffusions(element,physics);
            this.computeJacobians(element,physics);
            aux = (cell2mat(element.jacobians))*element.states(:);
            aux = reshape(aux,physics.equationCount,this.basisCount);
            element.residuals = aux - element.riemannR.*this.right' - element.riemannL.*this.left';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% WORKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            element.computeFluxesFromStates(physics);
            element.residuals = ...
                element.fluxes*this.gradientMatrix...
                - element.riemannR.*this.right'...
                - element.riemannL.*this.left';
            this.diffuseResiduals(element,physics);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            element.residuals = 2/element.dx * element.residuals ./ this.lumpedMassMatrixDiagonal;
            %%% element.residuals = 2/element.dx * element.residuals / this.massMatrix; % TESTING

        end
    end
    methods (Static)
        %% Low order residuals
        function diffuseResiduals(element,physics)
            % Applies AFC diffusion to the residual matrix, a la Möller &
            % Jaeschke, 2018 (eq. 16). Saves the diffusion block matrix
            % into the given element as a cell array.
            %
            % Precompute the "gradient differences" (Kuzmin et al. 2012, eq. 27):
            e = .5*abs(element.basis.gradientMatrix - element.basis.gradientMatrix');
            % Loop over edges:
            for edge = element.basis.edges
                % Aliases:
                r = edge(1); j = edge(2);
                % Compute a diffusion matrix block:
%%% ????????? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [D,L,R] = physics.getEigensystemAt(.5*(element.states(:,r) + element.states(:,j))); % Kuzmin et al. 2012, eq. 44
                D = e(r,j)*R*abs(D)*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.diffusions{r,j} = D;
                % Accumulate diffusion into the residuals:
                element.residuals(:,r) = element.residuals(:,r) +...
                    D*(element.states(:,j) - element.states(:,r));
            end
        end
        %% Initialize block matrices
        function initializeBlockMatrices(element,N,I)
            % Preallocates memory for the discrete Jacobian and artificial
            % diffusion block matrices on a given element, and computes the
            % connectivity matrix between blocks and block matrix entries.
            %
            % Blocks are arranged a la Kuzmin et al., 2012, defined so that
            % the following relationship between element residuals holds:
            % 
            %  R_Kuzmin == R(:)
            % 
            % where:
            %
            %  R_Kuzmin = element.jacobians*element.states(:)
            %  R = element.fluxes*this.gradientMatrix
            %
            % Arguments
            %  element: the patch to initialize
            %  N: number of basis functions in the patch
            %  I: number of system components
            %
            % Sparse cell array preallocation:
            element.jacobians = cell(N);
            element.jacobians(:) = {sparse(I,I)};
            element.diffusions = element.jacobians;
%%% DEPRECATED (fast, but too complicated to maintain %%%%%%%%%%%%%%%%%%%%%
%             % Connectivity matrix:
%             a = repmat(1:I,1,I);
%             b = repelem(1:I,I);
%             [r,j] = find(element.basis.massMatrix);
%             element.basis.blockConnectivity = (I*N*(I*(j-1)+b-1)+I*(r-1)+a)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %% Diffusion operator
        function computeDiffusions(element,physics)
            % Updates edge-by-edge the discrete numerical diffusion block 
            % matrix of an IGA patch.
            %
            % Precompute the "gradient differences" (Kuzmin et al. 2012, eq. 27):
            e = .5*abs(element.basis.gradientMatrix - element.basis.gradientMatrix');
            % Loop over edges (a la Kuzmin et al. 2012):
            for rj = element.basis.edges
                % Compute an "edge state":
                state = .5*(element.states(:,rj(1)) + element.states(:,rj(2))); % arithmetic (Kuzmin et al. 2012, eq. 44)
                % Evaluate flux jacobian eigensystem at the (r,j) edge:
                [D,L,R] = physics.getEigensystemAt(state);
                % Compute current diffusion block:
                D = e(rj(1),rj(2))*R*abs(D)*L; % Möller et al. 2018, eq. 16
                % Update numerical diffusion operator:
                element.diffusions{rj(1),rj(2)} = D;
                % Zero row sum condition for diagonal blocks (Kuzmin 2012, eq. 38):
                element.diffusions{rj(1),rj(1)} = element.diffusions{rj(1),rj(1)} - D;
            end
        end
        %% Jacobian operator
        function computeJacobians(element,physics)
            % Updates edge-by-edge the discrete Jacobian block matrix of 
            % an IGA patch. Applies existing diffusion matrix to it (low 
            % order Jacobian).
            %
            % Loop over pairs of control points with shared support:
            for rj = element.basis.modePairs
                % Evaluate flux jacobians at the j-th control point:
                K = physics.getJacobianAt(element.states(:,rj(2)));
                % Apply the discrete gradient operator:
                K = element.basis.gradientMatrix(rj(2),rj(1))*K;
                % Apply numerical diffusion:
                K = K + element.diffusions{rj(1),rj(2)};
                % Update low-order Jacobian operator:
                element.jacobians{rj(1),rj(2)} = K;
            end
        end
        %% Constrained L2 projection (vector)
        function project(mesh,limiter,fun,q)
            % Projects a function into a finite-dimensional approximation space
            % using AFC and FCT flux limiting (Zalesak's approach) to constrain
            % the approaximate solution to be L2-preserving and LED, while
            % maintaining a higher order than a lumped projection.
            %
            % DG coupling: candidate local extrema in non-interior modes
            % include the closest mode across the closest patch interface.
            %
            % Initialize mesh via unconstrained L2 projection:
            if nargin < 4
                Bspline.project(mesh,[],fun);
            else
                Bspline.project(mesh,[],fun,q);
            end
            
            % ...TO BE DONE...
%             arrayfun(@(element) DGIGA_AFC_vector.initializeBlockMatrices(...
%                 element,limiter.physics),mesh.elements); % preallocate jacobian and diffusion cell arrays
        
        end
        %% Constrained L2 projection (scalar)
        function projectScalar(mesh,~,fun,q)
            % Projects a function into a finite-dimensional approximation space
            % using AFC and FCT flux limiting (Zalesak's approach) to constrain
            % the approaximate solution to be L2-preserving and LED, while
            % maintaining a higher order than a lumped projection.
            %
            % DG coupling: candidate local extrema in non-interior modes
            % include the closest mode across the closest patch interface.
            %
            % Limitations (so far):
            %  - scalar physics only (nEqs = 1)
            %
            % Initialize mesh via unconstrained L2 projection:
            if nargin < 4
                Bspline.project(mesh,[],fun);
            else
                Bspline.project(mesh,[],fun,q);
            end
            % Preallocate memory for low order solutions within DG range:
            statesL = cell(1,3);
            statesL{2} = mesh.elements(1).states*mesh.elements(1).basis.massMatrix./mesh.elements(1).basis.lumpedMassMatrixDiagonal;
            % Apply AFC patch-wise, with inter-patch coupling:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                [nEqs,nB] = size(element.states);
                if nEqs ~= 1
                    error('Vector physics not supported. Yet.')
                end
                % Anti-diffusive fluxes (row: recieving mode; column: contributing mode):
                f = element.basis.massMatrix.*(element.states' - element.states);
                % Pre-limiting:
                f(f.*(statesL{2}-statesL{2}') > 0) = 0;
                % Net antidiffusive fluxes:
                Pp = full(sum(max(0,f),2));
                Pm = full(sum(min(-0,f),2));
                Pm(Pm == 0) = -0; % signed zeros to avoid -Inf later
                % Local extrema (within patch):
                extrema = repmat(statesL{2},nB,1);
                extrema(element.basis.massMatrix == 0) = nan;
                Qp = max(extrema,[],2);
                Qm = min(extrema,[],2);
                % Local extrema (across patch interfaces):
                if k > 1
                    i = 1:element.basis.degree; % left-most mode of right patch might be an extremum of the first p modes of this patch
                    Qp(i) = max(Qp(i),statesL{1}(end));
                    Qm(i) = min(Qm(i),statesL{1}(end));
                end
                if k < mesh.elementCount
                    statesL{3} = mesh.elements(k+1).states*mesh.elements(k+1).basis.massMatrix./mesh.elements(k+1).basis.lumpedMassMatrixDiagonal; % next patch's low order predictor
                    i = nB - (element.basis.degree-1:-1:0); % right-most mode of left patch might be an extremum of the last p modes of this patch
                    Qp(i) = max(Qp(i),statesL{3}(1));
                    Qm(i) = min(Qm(i),statesL{3}(1));
                end
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.limiterHistory.localMax = Qp;
                element.limiterHistory.localMin = Qm;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Distances to local extrema:
                Qp = element.basis.lumpedMassMatrixDiagonal'.*(Qp-statesL{2}');
                Qm = element.basis.lumpedMassMatrixDiagonal'.*(Qm-statesL{2}');
                % Modal correction factors:
                Rp = min(1,Qp./Pp);
                Rm = min(1,Qm./Pm);
                % Disable AFC at boundaries (Kuzmin et al, 2012; remark 5, pp. 163, bottom):
                if k == 1
                    Rp(1) = 1;
                    Rm(1) = 1;
                end
                if k == mesh.elementCount
                    Rp(end) = 1;
                    Rm(end) = 1;
                end
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.limiterHistory.Pp = Pp;
                element.limiterHistory.Pm = Pm;
                element.limiterHistory.Qp = Qp;
                element.limiterHistory.Qm = Qm;
                element.limiterHistory.Rp = Rp;
                element.limiterHistory.Rm = Rm;
                element.limiterHistory.Alpha = ones(nB,nB);
                [i,j] = find(f > 0); ids = sub2ind([nB nB],i,j);
                element.limiterHistory.Alpha(ids) = min(Rp(i),Rm(j));
                [i,j] = find(f < 0); ids = sub2ind([nB nB],i,j);
                element.limiterHistory.Alpha(ids) = min(Rm(i),Rp(j));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Limit the antidiffusive fluxes:
                [i,j] = find(f > 0); ids = sub2ind([nB nB],i,j); % indices of positive fluxes
                f(ids) = min(Rp(i),Rm(j)).*f(ids);
                [i,j] = find(f < 0); ids = sub2ind([nB nB],i,j); % indices of negative fluxes
                f(ids) = min(Rm(i),Rp(j)).*f(ids);
                % Use low-order predictor and limited antidiffusive fluxes:
                element.states = statesL{2} + sum(f,2)'./element.basis.lumpedMassMatrixDiagonal;
                % Cycle low order solutions along their cell array:
                statesL(1:2) = statesL(2:3);
            end
        end
        function project_Matthias(mesh,~,fun,q)
            % Projects a function into a finite-dimensional approximation space
            % using AFC and FCT flux limiting (Zalesak's approach) to constrain
            % the approaximate solution to be L2-preserving and LED, while
            % maintaining a higher order than a lumped projection.
            %
            % Does the DG coupling such that P and Q are unique at patch
            % interfaces (and, therefore, so is R).
            %
            % Limitations (so far):
            %  - scalar physics only (nEqs = 1)
            %
            % Initialize mesh via unconstrained L2 projection:
            if nargin < 4
                Bspline.project(mesh,[],fun);
            else
                Bspline.project(mesh,[],fun,q);
            end
            % Preallocate some memory:
            f = cell(1,mesh.elementCount);
            statesL = cell(1,mesh.elementCount);
            Pp = cell(1,mesh.elementCount);
            Pm = cell(1,mesh.elementCount);
            Qp = cell(1,mesh.elementCount);
            Qm = cell(1,mesh.elementCount);
            % Sart Zalesak's algorithm patch-wise:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                [nEqs,nB] = size(element.states);
                if nEqs ~= 1
                    error('Vector physics not supported. Yet.')
                end
                % Anti-diffusive fluxes (row: recieving mode; column: contributing mode):
                f{k} = element.basis.massMatrix.*(element.states' - element.states);
                % Pre-limiting:
                statesL{k} = element.states*element.basis.massMatrix./element.basis.lumpedMassMatrixDiagonal; % low order predictor
                f{k}(f{k}.*(statesL{k}-statesL{k}') > 0) = 0;
                % Net antidiffusive fluxes:
                Pp{k} = full(sum(max(0,f{k}),2));
                Pm{k} = full(sum(min(-0,f{k}),2));
                Pm{k}(Pm{k} == 0) = -0; % signed zeros to avoid -Inf later
                % Local extrema (within patch):
                extrema = repmat(statesL{k},nB,1);
                extrema(element.basis.massMatrix == 0) = nan;
                % Distances to local extrema:
                Qp{k} = element.basis.lumpedMassMatrixDiagonal'.*(max(extrema,[],2)-statesL{k}');
                Qm{k} = element.basis.lumpedMassMatrixDiagonal'.*(min(extrema,[],2)-statesL{k}');
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             element.limiterHistory.localMax = max(extrema,[],2);
                %             element.limiterHistory.localMin = min(extrema,[],2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            % Perform inter-patch coupling of P and Q variables:
            k = 1; % skip 1st patch
            for element = mesh.elements(2:end)
                k = k + 1;
                % P+-:
                Pp{k}(1) = max(Pp{k-1}(end),Pp{k}(1));
                Pp{k-1}(end) = Pp{k}(1); % left interface of patch k <-> right of patch k-1
                Pm{k}(1) = min(Pm{k-1}(end),Pm{k}(1));
                Pm{k-1}(end) = Pm{k}(1);
                % Q+-:
                Qp{k}(1) = max(Qp{k-1}(end),Qp{k}(1));
                Qp{k-1}(end) = Qp{k}(1);
                Qm{k}(1) = min(Qm{k-1}(end),Qm{k}(1));
                Qm{k-1}(end) = Qm{k}(1);
            end
            % Finish with Zalesak patch-wise:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                % Modal correction factors:
                Rp = min(1,Qp{k}./Pp{k});
                Rm = min(1,Qm{k}./Pm{k});
                % Disable AFC at boundaries (Kuzmin et al, 2012; remark 5, pp. 163, bottom):
                if k == 1
                    Rp(1) = 1;
                    Rm(1) = 1;
                end
                if k == mesh.elementCount
                    Rp(end) = 1;
                    Rm(end) = 1;
                end
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             element.limiterHistory.Pp = Pp{k};
                %             element.limiterHistory.Pm = Pm{k};
                %             element.limiterHistory.Qp = Qp{k};
                %             element.limiterHistory.Qm = Qm{k};
                %             element.limiterHistory.Rp = Rp;
                %             element.limiterHistory.Rm = Rm;
                %             element.limiterHistory.Alpha = ones(nB,nB);
                %             [i,j] = find(f{k} > 0); ids = sub2ind([nB nB],i,j);
                %             element.limiterHistory.Alpha(ids) = min(Rp(i),Rm(j));
                %             [i,j] = find(f{k} < 0); ids = sub2ind([nB nB],i,j);
                %             element.limiterHistory.Alpha(ids) = min(Rm(i),Rp(j));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Limit the antidiffusive fluxes:
                [i,j] = find(f{k} > 0); ids = sub2ind([nB nB],i,j); % indices of positive fluxes
                f{k}(ids) = min(Rp(i),Rm(j)).*f{k}(ids);
                [i,j] = find(f{k} < 0); ids = sub2ind([nB nB],i,j); % indices of negative fluxes
                f{k}(ids) = min(Rm(i),Rp(j)).*f{k}(ids);
                % Use low-order predictor and limited antidiffusive fluxes:
                element.states = statesL{k} + sum(f{k},2)'./element.basis.lumpedMassMatrixDiagonal;
            end
        end
    end
end