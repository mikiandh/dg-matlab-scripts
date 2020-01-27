classdef Wang_vector < Limiter
    properties
        everyStage = true;
        physics
        beta2 = 0.75;
    end
    methods
        %% Constructor
        function this = Wang_vector(physics)
            this.physics = physics;
        end
        %% AP-TVD detector + PFGM limiter (Wang 2009) - partially vectorized -
        function apply(this,mesh,~)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some quantities:
            eqnIds = 1:this.physics.equationCount; % precomputed array of state vector component ids
            cellIds = 1:mesh.elementCount;
            q = mesh.getNodalCoeffs(eqnIds);
            Q = mesh.getModalCoeffs(eqnIds);
            [L,R] = this.physics.getEigenvectorsAt(Q(:,1,:));
            Q_hat = Q;
            for k = cellIds
                mesh.elements(k).limiterHistory = zeros(length(eqnIds),1);
                Q_hat(:,:,k) = L(:,:,k)*Q_hat(:,:,k);
            end
            % AP-TVD marker:
            TROUBLED_COMPONENTS = false(length(eqnIds),length(cellIds)); % troubled state vector components in each element
            % Step 1:
            qMax = 1.001*max([Q(:,1,1:end-2) Q(:,1,2:end-1) Q(:,1,3:end)],[],2);
            qMin = 0.999*min([Q(:,1,1:end-2) Q(:,1,2:end-1) Q(:,1,3:end)],[],2);
            for k = cellIds(2:end-1)
                TROUBLED_COMPONENTS(:,k) = any(q(:,:,k) > qMax(:,:,k-1) | q(:,:,k) < qMin(:,:,k-1),2);
            end
            TROUBLED_CELLS = any(TROUBLED_COMPONENTS,1);
            % Step 2:
            for k = cellIds(TROUBLED_CELLS)
                diffL = this.beta2*L(:,:,k)*(Q(:,1,k+1) - Q(:,1,k));
                diffR = this.beta2*L(:,:,k)*(Q(:,1,k) - Q(:,1,k-1));
                TROUBLED_COMPONENTS(:,k) = this.minmod1D(Q_hat(:,2,k),diffR,diffL) ~= Q_hat(:,2,k);
            end
            TROUBLED_CELLS = any(TROUBLED_COMPONENTS,1);
            if all(~TROUBLED_CELLS)
                return
            end
            TROUBLED_COMPONENTS = TROUBLED_COMPONENTS';
            % PFGM limiter:
            for j = mesh.maxDegree:-1:1
                aux = this.beta2/(2*j-1);
                for k = cellIds(TROUBLED_CELLS)
                    for i = eqnIds(TROUBLED_COMPONENTS(k,:))
                        % Characteristic decomposition (type 2):
                        moment = Q_hat(i,j+1,k);
                        diffR = aux*L(i,:,k)*(Q(:,j,k+1) - Q(:,j,k));
                        diffL = aux*L(i,:,k)*(Q(:,j,k) - Q(:,j,k-1));
                        Q_hat(i,j+1,k) = this.minmod(moment,diffR,diffL);
                        if Q_hat(i,j+1,k) == moment && moment ~= 0
                            TROUBLED_COMPONENTS(k,i) = false; % component is no longer troubled
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