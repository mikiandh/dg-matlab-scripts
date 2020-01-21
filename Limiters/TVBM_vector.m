classdef TVBM_vector < Limiter
    properties
        everyStage = true;
        M % modified Minmod function parameter
        physics
    end
    methods
        %% Constructor
        function this = TVBM_vector(physics,M)
            if nargin == 0
                M = 0; % TVDM
            end
            this.M = M;
            this.physics = physics;
        end
        %% TVBM limiter (simple local characteristic decomposition)
        function apply(this,mesh)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some quantities:
            eqnIds = 1:this.physics.equationCount; % precomputed array of state vector component ids
            cellIds = 1:mesh.elementCount; % precomputed array of cell ids
            Q = mesh.getModalCoeffs(eqnIds);
            [L,R] = this.physics.getEigenvectorsAt(Q(:,1,:));
            Q_hat = Q;
            for k = cellIds
                mesh.elements(k).limiterHistory = zeros(length(eqnIds),1);
                Q_hat(:,:,k) = L(:,:,k)*Q_hat(:,:,k);
            end
            % Do the limiting:
            for i = eqnIds
                for k = cellIds(2:end-1)
                    tilda = Q_hat(i,2,k);
                    diffR = L(i,:,k)*(Q(:,1,k+1) - Q(:,1,k));
                    diffL = L(i,:,k)*(Q(:,1,k) - Q(:,1,k-1));
                    m = this.M*mesh.elements(k).dx^2;
                    Q_hat(i,2,k) = this.modminmod(tilda,diffR,diffL,m);
                    if Q_hat(i,2,k) ~= tilda
                        Q_hat(i,3:end,k) = 0;
                        mesh.elements(k).limiterHistory(i) = mesh.elements(k).basis.degree;
                    end
                end
            end
            % Revert characteristic decomposition:
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
            % Overwrite solution with the limited coefficients:
            mesh.setModalCoeffs(Q,eqnIds);
        end
    end
    methods (Static)
        %% Modified minmod function (3 scalar arguments)
        function d = modminmod(a,b,c,m)
            d = sign(a);
            if abs(a) <= m
                d = a;
            elseif d == sign(b) && d == sign(c)
                d = d*min(abs([a b c]));
            else
                d = 0;
            end
        end
    end
end