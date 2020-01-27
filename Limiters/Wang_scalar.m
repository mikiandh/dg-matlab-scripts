classdef Wang_scalar < Limiter
    properties
        everyStage = true;
        beta2 = 0.75;
    end
    methods
        %% AP-TVD detector + PFGM limiter (Wang 2009)
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
            % AP-TVD marker:
            TROUBLED = false(size(cellIds));
            for k = cellIds(2:end-1)
                means = [Q(1,1,k-1) Q(1,1,k) Q(1,1,k+1)];
                qMax = max(means);
                qMin = min(means);
                mesh.elements(k).interpolateStateAtEdges;
                q = [mesh.elements(k).stateL mesh.elements(k).states mesh.elements(k).stateR];
                if any(q > 1.001*qMax | q < 0.999*qMin)
                    diffs = this.beta2*diff(means);
                    if this.minmod(Q(1,2,k),diffs(2),diffs(1)) ~= Q(1,2,k)
                        TROUBLED(k) = true; % current cell is troubled
                    end
                end
            end
            % PFGM limiter:
            for j = mesh.maxDegree:-1:1 % basis component
                aux = this.beta2/(2*j-1);
                for k = cellIds(TROUBLED) % troubled element
                    moment = Q(1,j+1,k);
                    diffR = Q(1,j,k+1) - Q(1,j,k);
                    diffL = Q(1,j,k) - Q(1,j,k-1);
                    Q(1,j+1,k) = this.minmod(moment,aux*diffR,aux*diffL);
                    if Q(1,j+1,k) == moment && moment ~= 0
                        TROUBLED(k) = 0; % cell is no longer troubled
                    else
                        mesh.elements(k).limiterHistory = mesh.elements(k).limiterHistory + 1;
                    end
                end
            end
            % Overwrite solution with the limited modal coefficients:
            mesh.setModalCoeffs(Q,1);
        end
    end
end