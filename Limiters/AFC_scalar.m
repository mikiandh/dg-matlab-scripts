classdef AFC_scalar < Limiter
    properties
        everyStage = false; % false -> too diffusive; true -> too antidiffusive
        physics
    end
    methods
        %% Constructor
        function this = AFC_scalar(physics)
            this.physics = physics;
        end
        %% Algebraic Flux Correction (no interpatch coupling)
        function applyNONE(this,mesh,timeDelta)
            % Employs the AFC procedure to apply FCT limiting to the modes
            % of a (quasi-)nodal finite-element discretization.
            %
            % Assumes that the mesh contains the low order solution and 
            % residual (i.e. time-derivative) at the next time-step
            % ("transported and diffused"). It then reconstructs a limited 
            % high order solution (which is LED) via an AFC-based FCT 
            % procedure applied to the antidiffusive flux components 
            % associated to each combination of two modes.
            %
            % Assumes a scalar conservation law.
            %
            % Safety check:
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Compute the low-order modal time-derivatives:
            mesh.computeResiduals(this.physics);
            % Apply AFC patch-wise, with inter-patch coupling:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                [nEqs,nB] = size(element.states);
                if nEqs ~= 1
                    error('Vector physics not supported. Yet.')
                end
                % Anti-diffusive fluxes (row: recieving mode; column: contributing mode):
                diffusionMatrix = this.physics.getConvectionMatrix(element.states,element.basis);
                diffusionMatrix = element.basis.applyDiffusion(diffusionMatrix) - diffusionMatrix;
                f = element.basis.massMatrix.*(element.residuals' - element.residuals) + diffusionMatrix.*(element.states' - element.states);
                % Pre-limiting:
                f(f.*(element.states-element.states') > 0) = 0;
                % Net antidiffusive fluxes:
                Pp = full(sum(max(0,f),2));
                Pm = full(sum(min(-0,f),2));
                Pm(Pm == 0) = -0; % signed zeros to avoid -Inf later
                % Local extrema (within patch):
                extrema = repmat(element.states,nB,1);
                extrema(element.basis.massMatrix == 0) = nan; % "mask-out" the modes outside local influence
                Qp = max(extrema,[],2);
                Qm = min(extrema,[],2);
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.limiterHistory.localMax = Qp;
                element.limiterHistory.localMin = Qm;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Distances to local extrema:
                Qp = element.basis.lumpedMassMatrixDiagonal'.*(Qp-element.states')/timeDelta;
                Qm = element.basis.lumpedMassMatrixDiagonal'.*(Qm-element.states')/timeDelta;
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
                % Explicit correction using limited antidiffusive fluxes:
                element.states = element.states + timeDelta*sum(f,2)'./element.basis.lumpedMassMatrixDiagonal;
            end
        end
        %% Algebraic Flux Correction
        function applyOLD(this,mesh,timeDelta)
            % Employs the AFC procedure to apply FCT limiting to the modes
            % of a (quasi-)nodal finite-element discretization.
            %
            % Assumes that the mesh contains the low order solution and 
            % residual (i.e. time-derivative) at the next time-step
            % ("transported and diffused"). It then reconstructs a limited 
            % high order solution (which is LED) via an AFC-based FCT 
            % procedure applied to the antidiffusive flux components 
            % associated to each combination of two modes.
            %
            % Assumes a scalar conservation law.
            %
            % Safety check:
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Compute the low-order modal time-derivatives:
            mesh.computeResiduals(this.physics);
            % Apply AFC patch-wise, with inter-patch coupling:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                [nEqs,nB] = size(element.states);
                if nEqs ~= 1
                    error('Vector physics not supported. Yet.')
                end
                % Anti-diffusive fluxes (row: recieving mode; column: contributing mode):
                diffusionMatrix = this.physics.getConvectionMatrix(element.states,element.basis);
                diffusionMatrix = element.basis.applyDiffusion(diffusionMatrix) - diffusionMatrix;
                f = element.basis.massMatrix.*(element.residuals' - element.residuals) + diffusionMatrix.*(element.states' - element.states);
                % Pre-limiting:
                f(f.*(element.states-element.states') > 0) = 0;
                % Net antidiffusive fluxes:
                Pp = full(sum(max(0,f),2));
                Pm = full(sum(min(-0,f),2));
                Pm(Pm == 0) = -0; % signed zeros to avoid -Inf later
                % Local extrema (within patch):
                extrema = repmat(element.states,nB,1);
                extrema(element.basis.massMatrix == 0) = nan; % "mask-out" the modes outside local influence
                Qp = max(extrema,[],2);
                Qm = min(extrema,[],2);
                % Local extrema (across patch interfaces):
                if k > 1
                    i = 1:element.basis.degree; % left-most mode of right patch might be an extremum of the first p modes of this patch
                    Qp(i) = max(Qp(i),stateLowLeftRight);
                    Qm(i) = min(Qm(i),stateLowLeftRight);
                end
                if k < mesh.elementCount
                    i = nB - (element.basis.degree-1:-1:0); % right-most mode of left patch might be an extremum of the last p modes of this patch
                    Qp(i) = max(Qp(i),mesh.elements(k+1).states(1));
                    Qm(i) = min(Qm(i),mesh.elements(k+1).states(1));
                end
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.limiterHistory.localMax = Qp;
                element.limiterHistory.localMin = Qm;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Distances to local extrema:
                Qp = element.basis.lumpedMassMatrixDiagonal'.*(Qp-element.states')/timeDelta;
                Qm = element.basis.lumpedMassMatrixDiagonal'.*(Qm-element.states')/timeDelta;
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
                % Store the low-order coefficient of this patch's rightmost mode:
                stateLowLeftRight = element.states(end); % low-order, patch on the left of upcoming patch, righ-most mode
                % Explicit correction using limited antidiffusive fluxes:
                element.states = element.states + timeDelta*sum(f,2)'./element.basis.lumpedMassMatrixDiagonal;
            end
        end
        %% Algebraic Flux Correction (P & Q coupling version)
        function apply(this,mesh,timeDelta)
            % Employs the AFC procedure to apply FCT limiting to the modes
            % of a (quasi-)nodal finite-element discretization. Interpatch
            % coupling is made via P and Q communication across edges.
            %
            % Assumes that the mesh contains the low order solution and 
            % residual (i.e. time-derivative) at the next time-step
            % ("transported and diffused"). It then reconstructs a limited 
            % high order solution (which is LED) via an AFC-based FCT 
            % procedure applied to the antidiffusive flux components 
            % associated to each combination of two modes.
            %
            % Assumes a scalar conservation law.
            %
            % Safety check:
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Compute the low-order modal time-derivatives:
            mesh.computeResiduals(this.physics);
            % Preallocations:
            f = cell(1,mesh.elementCount);
            Pp = cell(1,mesh.elementCount);
            Pm = cell(1,mesh.elementCount);
            Qp = cell(1,mesh.elementCount);
            Qm = cell(1,mesh.elementCount);
            % Apply AFC patch-wise, with inter-patch coupling:
            k = 0;
            for element = mesh.elements
                k = k + 1;
                [nEqs,nB] = size(element.states);
                if nEqs ~= 1
                    error('Vector physics not supported. Yet.')
                end
                % Anti-diffusive fluxes (row: recieving mode; column: contributing mode):
                diffusionMatrix = this.physics.getConvectionMatrix(element.states,element.basis);
                diffusionMatrix = element.basis.applyDiffusion(diffusionMatrix) - diffusionMatrix;
                f{k} = element.basis.massMatrix.*(element.residuals' - element.residuals) + diffusionMatrix.*(element.states' - element.states);
                % Pre-limiting:
                f{k}(f{k}.*(element.states-element.states') > 0) = 0;
                % Net antidiffusive fluxes:
                Pp{k} = full(sum(max(0,f{k}),2));
                Pm{k} = full(sum(min(-0,f{k}),2));
                Pm{k}(Pm{k} == 0) = -0; % signed zeros to avoid -Inf later
                % Local extrema (within patch):
                extrema = repmat(element.states,nB,1);
                extrema(element.basis.massMatrix == 0) = nan; % "mask-out" the modes outside local influence
                Qp{k} = max(extrema,[],2);
                Qm{k} = min(extrema,[],2);
                % Distances to local extrema (within patch):
                Qp{k} = element.basis.lumpedMassMatrixDiagonal'.*(Qp{k}-element.states')/timeDelta;
                Qm{k} = element.basis.lumpedMassMatrixDiagonal'.*(Qm{k}-element.states')/timeDelta;
                
                %%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                element.limiterHistory.localMax = max(extrema,[],2);
                element.limiterHistory.localMin = min(extrema,[],2);
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
                [~,nB] = size(element.states);
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
                element.limiterHistory.Pp = Pp{k};
                element.limiterHistory.Pm = Pm{k};
                element.limiterHistory.Qp = Qp{k};
                element.limiterHistory.Qm = Qm{k};
                element.limiterHistory.Rp = Rp;
                element.limiterHistory.Rm = Rm;
                element.limiterHistory.Alpha = ones(nB,nB);
                [i,j] = find(f{k} > 0); ids = sub2ind([nB nB],i,j);
                element.limiterHistory.Alpha(ids) = min(Rp(i),Rm(j));
                [i,j] = find(f{k} < 0); ids = sub2ind([nB nB],i,j);
                element.limiterHistory.Alpha(ids) = min(Rm(i),Rp(j));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Limit the antidiffusive fluxes:
                [i,j] = find(f{k} > 0); ids = sub2ind([nB nB],i,j); % indices of positive fluxes
                f{k}(ids) = min(Rp(i),Rm(j)).*f{k}(ids);
                [i,j] = find(f{k} < 0); ids = sub2ind([nB nB],i,j); % indices of negative fluxes
                f{k}(ids) = min(Rm(i),Rp(j)).*f{k}(ids);
                % Explicit correction using limited antidiffusive fluxes:
                element.states = element.states + timeDelta*sum(f{k},2)'./element.basis.lumpedMassMatrixDiagonal;
            end
        end
    end
end