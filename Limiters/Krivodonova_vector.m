classdef Krivodonova_vector < Limiter
    properties
        everyStage = true;
        physics
        ratio
    end
    methods
        %% Constructor
        function this = Krivodonova_vector(physics,ratio)
            this.physics = physics;
            this.ratio = ratio;
        end
        %% High order moment limiter (Krivodonova 2007)
        function apply(this,mesh,~)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some things:
            eqnIds = 1:this.physics.equationCount; % precomputed array of state vector component ids
            cellIds = 1:mesh.elementCount; % precomputed array of cell ids
            Q = mesh.getModalCoeffs(eqnIds);
            [L,R] = this.physics.getEigenvectorsAt(Q(:,1,:));
            Q_hat = Q;
            for k = cellIds
                mesh.elements(k).limiterHistory = zeros(length(eqnIds),1); % reset limiter history
                Q_hat(:,:,k) = L(:,:,k)*Q_hat(:,:,k);
            end
            % Do the limiting:
            TROUBLED_COMPONENTS = true(length(cellIds),length(eqnIds)); % troubled state vector components in each element
            TROUBLED_CELLS = true(1,length(cellIds));
            TROUBLED_CELLS(:,[1 end]) = 0; % exclude elements adjacent to boundaries
            for j = mesh.maxDegree:-1:1
                aux = this.ratio*(0.5/(2*j-1) - 1) + 1;
                for k = cellIds(TROUBLED_CELLS)
                    for i = eqnIds(TROUBLED_COMPONENTS(k,:))
                        % Characteristic decomposition (type 2):
                        moment = Q_hat(i,j+1,k);
                        diffR = aux*L(i,:,k)*(Q(:,j,k+1) - Q(:,j,k));
                        diffL = aux*L(i,:,k)*(Q(:,j,k) - Q(:,j,k-1));
                        Q_hat(i,j+1,k) = this.minmod(moment,diffR,diffL);
                        if Q_hat(i,j+1,k) == moment && moment ~= 0
                            TROUBLED_COMPONENTS(k,i) = 0; % component is no longer troubled
                        else
                            mesh.elements(k).limiterHistory(i) = mesh.elements(k).limiterHistory(i) + 1;
                        end
                    end
                    TROUBLED_CELLS(k) = any(TROUBLED_COMPONENTS(k,:)); % cell is troubled if any of its state vector components is troubled
                end
            end
            % Transform limited modal coefficients back to conserved variables:
            for k = cellIds
                Q(:,:,k) = R(:,:,k)*Q_hat(:,:,k);
            end
            % Limit down to p = 0 where necessary:
            [legendreL,legendreR] = Legendre.getComponentsAtEdges(mesh.maxDegree);
            for k = cellIds(2:end-1)
                state1 = Q(:,:,k-1)*legendreR;
                state2 = Q(:,:,k)*legendreL;
                state3 = Q(:,:,k)*legendreR;
                state4 = Q(:,:,k+1)*legendreL;
                if ~this.physics.checkWithinPhysicalBounds(state1,state2) || ~this.physics.checkWithinPhysicalBounds(state3,state4)
                    Q_hat(:,2:end,k) = 0;
                    Q(:,:,k) = R(:,:,k)*Q_hat(:,:,k);
                end
            end
            % Overwrite solution with the limited modal coefficients:
            mesh.setModalCoeffs(Q,eqnIds);
        end
    end
end