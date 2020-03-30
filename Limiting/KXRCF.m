classdef KXRCF < Sensor
    %
    % Sensor by Krivodonova et al., 2004. Based on the superconvergence 
    % property that DG has at outflow boundaries of cells where the exact
    % solution is smooth, where order of accuracy rises to 2(p+1).
    %
    % Applied characteristic-wise.
    %
    methods
        %% Sensor
        function apply(this,mesh,solver)
            % Apply default sensor first:
            apply@Sensor(this,mesh);
            % Evaluate state at all edges:
            for element = mesh.elements
                element.interpolateStateAtEdges;
            end
            % Check all possibly troubled elements:
            for element = mesh.elements([mesh.elements.isTroubled])
                % Precompute some stuff:
                ref = abs(element.basis.getLegendre(element,1))*element.dx^(.5*element.dofCount);
                qL = element.stateL; % left edge of element k
                qLL = element.elementL.stateR; % right edge of element k-1
                qR = element.stateR; % right edge of element k
                qRR = element.elementR.stateL; % left edge of element k+1
                % Test left edge:
                try
                    [D,L] = solver.physics.getEigensystemAt(qLL,qL);
                catch
                    continue
                end
                i = diag(D) < 0; % out-going characteristics
                if any(L(i,:)*abs(qLL - qL)./(L(i,:)*ref) > 1)
                    continue
                end
                % Test right edge:
                try
                    [D,L] = solver.physics.getEigensystemAt(qR,qRR);
                catch
                    continue
                end
                i = diag(D) > 0; % out-going characteristics
                if any(L(i,:)*abs(qRR - qR)./(L(i,:)*ref) > 1)
                    continue
                end
                % If it got this far, element is untroubled:
                element.isTroubled = false;
            end
        end
    end
end