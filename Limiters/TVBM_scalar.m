classdef TVBM_scalar < Limiter
    properties
        everyStage = true;
        M % modified Minmod function parameter
    end
    methods
        %% Constructor
        function this = TVBM_scalar(M)
            if nargin == 0
                M = 0; % TVDM
            end
            this.M = M;
        end
        %% TVBM limiter
        function apply(this,mesh,~)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some quantities:
            cellIds = 1:mesh.elementCount;
            for element = mesh.elements
                element.interpolateStateAtEdges;
                element.limiterHistory = 0;
            end
            % Extract modal coefficients:
            Q = mesh.getModalCoeffs(1);
            % Do the limiting:
            for k = cellIds(2:end-1)
                tilda = Q(1,2,k);
                diffR = Q(1,1,k+1) - Q(1,1,k);
                diffL = Q(1,1,k) - Q(1,1,k-1);
                m = this.M*mesh.elements(k).dx^2;
                Q(1,2,k) = this.modminmod(tilda,diffR,diffL,m);
                if Q(1,2,k) ~= tilda
                    Q(1,3:end,k) = 0;
                    mesh.elements(k).limiterHistory = mesh.maxDegree;
                end
            end
            % Overwrite solution with the limited coefficients:
            mesh.setModalCoeffs(Q,1);
        end
    end
    methods (Static)
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