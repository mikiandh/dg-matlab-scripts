classdef Krivodonova_scalar < Limiter
    properties
        everyStage = true;
        ratio
    end
    methods
        %% Constructor
        function this = Krivodonova_scalar(ratio)
            this.ratio = ratio;
        end
        %% High order moment limiter (Krivodonova 2007)
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
            % Hierarchical moment limiting:
            TROUBLED = true(size(cellIds));
            TROUBLED([1 end]) = 0; % exclude elements adjacent to boundaries
            for j = mesh.maxDegree:-1:1 % basis component
                aux = this.ratio*(0.5/(2*j-1) - 1) + 1;
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