classdef Burbeau_scalar < Limiter
    properties
        everyStage = true;
    end
    methods
        %% Burbeau et al. moment limiter (2001)
        function apply(this,mesh,~)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some quantities:
            cellIds = 1:mesh.elementCount;
            Q = mesh.getModalCoeffs(1);
            for k = cellIds
                mesh.elements(k).limiterHistory = 0;
            end
            % Do the limiting:
            TROUBLED = true(size(cellIds));
            TROUBLED([1 end]) = 0; % exclude elements adjacent to boundaries
            Q_mod = Q; % limited moments
            for j = mesh.maxDegree:-1:1
                aux1 = 2*j-1;
                aux2 = 1/aux1;
                for k = cellIds(TROUBLED)
                    moment = aux1*Q(1,j+1,k);
                    diffR = Q(1,j,k+1) - Q(1,j,k);
                    diffL = Q(1,j,k) - Q(1,j,k-1);
                    Q_m = aux2*this.minmod(moment,diffR,diffL);
                    if abs(Q(1,j+1,k) - Q_m) < 2*eps && moment ~= 0 % if Q_m == Q(1,j+1,k) && moment ~= 0
                        TROUBLED(k) = 0; % cell is no longer troubled
                        continue
                    end
                    mesh.elements(k).limiterHistory = mesh.elements(k).limiterHistory + 1;
                    diffR = Q(1,j,k+1) - aux1*Q(1,j+1,k+1);
                    diffR = diffR - Q(1,j,k);
                    diffL = Q(1,j,k-1) - aux1*Q(1,j+1,k-1);
                    diffL = Q(1,j,k) - diffL;
                    Q_mod(1,j+1,k) = aux2*this.minmod(moment,diffR,diffL);
                    Q_mod(1,j+1,k) = this.maxmod(Q_m,Q_mod(1,j+1,k));
                end
            end
            % Overwrite solution with the limited modal coefficients:
            mesh.setModalCoeffs(Q_mod,1);
        end
    end
    methods (Static)
        function c = maxmod(a,b)
            % maxmod function for 2 scalar arguments (Burbeau et al 2001)
            s = sign(a);
            if s == sign(b)
                c = s*max(abs([a b]));
            else
                c = 0;
            end
        end
    end
end