classdef Biswas_scalar < Limiter
    properties
        everyStage = true;
    end
    methods
        %% Biswas et al. moment limiter (1994)
        function apply(this,mesh,~)
            if mesh.maxDegree == 0
                return % nothing to limit
            end
            % Precompute some quantities:
            cellIds = 1:mesh.elementCount;
            Q = mesh.getModalCoeffs(1);
            %%% DEBUGGING %%%
%             disp 'Modes before limiting:'
%             i = 1:mesh.maxDegree+1;
%             disp(reshape(Q(1,i,:),length(i),[]))
%             disp(Q(1,2,8))
            %%%%%%%%%%%%%%%%%
            for k = cellIds
                mesh.elements(k).limiterHistory = 0;
            end
            % Hierarchical moment limiting:
            for l = 1:2 % sweeps (Biswas suggests 2)
                TROUBLED = true(size(cellIds));
                TROUBLED([1 end]) = 0; % exclude elements adjacent to boundaries
                for j = mesh.maxDegree:-1:1 % basis component
                    aux1 = 2*j-1;
                    aux2 = 1/aux1;
                    for k = cellIds(TROUBLED) % troubled element
                        moment = Q(1,j+1,k);
                        diffR = Q(1,j,k+1) - Q(1,j,k);
                        diffL = Q(1,j,k) - Q(1,j,k-1);
                        if 0
                            disp('moment:')
                            disp(moment)
                            disp('diffR:')
                            disp(diffR)
                            disp('diffL:')
                            disp(diffL)
                        end
                        Q(1,j+1,k) = aux2*this.minmod(aux1*moment,diffR,diffL);
                        if abs(Q(1,j+1,k) - moment) < 2*eps && moment ~= 0 % multiplication and division by 'aux' introduces round-off errors
                            TROUBLED(k) = 0; % cell is no longer troubled
                        else
                            mesh.elements(k).limiterHistory = mesh.elements(k).limiterHistory + 1;
                        end
                    end
                end
            end
            % Overwrite solution with the limited modal coefficients:
            mesh.setModalCoeffs(Q,1);
            %%% DEBUGGING %%%
%             disp 'Modes after limiting:'
%             disp(reshape(Q(1,i,:),length(i),[]))
%             disp(Q(1,2,8))
            %%%%%%%%%%%%%%%%%
        end
    end
end