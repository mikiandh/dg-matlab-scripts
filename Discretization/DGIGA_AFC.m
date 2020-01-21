classdef DGIGA_AFC < Bspline
    properties
        lumpedMassMatrixDiagonal % non-zero entries of the lumped mass matrix
    end
    methods
        %% Constructor
        function this = DGIGA_AFC(varargin)
            this@Bspline(varargin{:});
            if nargin > 0
                this.lumpedMassMatrixDiagonal = full(sum(this.massMatrix,1));
            end
        end
        %% Instantiate from prototype
        function this = clone(prototype,degree)
            this = DGIGA_AFC(prototype.knots,degree,prototype.smoothness);
        end
        %% Assemble convection operator
        function convectionMatrix = getConvectionMatrix(this,states,physics)
            % Assembles the high order convection matrix based on equation 
            % 23 of Kuzmin et al. (2012).
            %
%%% NO IDEA HOW TO GENERALIZE THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A = physics.getJacobian(states);
            convectionMatrix = this.gradientMatrix.*A(:); % convection matrix -> K_ij = A_j � c_ji
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %% DGIGA-AFC (low order) operator
        function computeResiduals(this,element,physics)
            element.computeFluxesFromStates(physics);
            convection = this.getConvectionMatrix(element.states,physics);
            convection = this.applyDiffusion(convection);
            element.residuals = 2/element.dx*(...
                element.states*convection - ...
                element.riemannR.*this.right' - ...
                element.riemannL.*this.left');
            element.residuals = element.residuals ./ this.lumpedMassMatrixDiagonal;
        end
    end
    methods (Static)
        %% Diffuse the convection operator
        function convectionMatrix = applyDiffusion(convectionMatrix)
            % Applies artificial diffusion to the discrete convection 
            % operator.
            %
            diffusion = full(max(-convectionMatrix,-convectionMatrix')); %%% OUT OF IDEAS
            diffusion = max(diffusion,0);
            mask = logical(eye(size(diffusion)));
            % Discrete upwinding (i.e. adding artificial diffusion):
            diffusion(mask) = -sum(diffusion,2) + diffusion(mask);
            convectionMatrix = convectionMatrix + diffusion;
            convectionMatrix = sparse(convectionMatrix);                 %%%
        end
        %% Constrained L2 projection
        function project(mesh,~,fun,q)
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