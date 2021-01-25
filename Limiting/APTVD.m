classdef APTVD < Sensor
    %
    % Accuracy-preserving TVD sensor based on Wang, 2009. Identifies local 
    % smoothness in the first derivatives of the solution (indicating a
    % smooth extrema). Implemented using Legendre coefficients in lieu of
    % cell-averaged derivatives (motivated by Krivodonova 2007).
    %
    methods
        %% Sensor
        function apply(this,mesh,solver,priority)
            % Apply default sensor first:
            apply@Sensor(this,mesh,solver,priority);
            % Mark all cells containing one or more local extrema:
            isTroubled = [mesh.elements.isTroubled];
            for element = mesh.elements(isTroubled(:,:,priority))
                % Get min/max cell average array:
                avgs = [element.getLegendre(1)...
                    element.elementL.getLegendre(1)...
                    element.elementR.getLegendre(1)];
                % Get solution samples:
                samples = element.interpolateStateAtCoords([element.xL element.getDofCoords' element.xR]);
                % Unmark if no local extremum is present:
                element.isTroubled(:,:,priority) = any(any(samples > 1.001*max(avgs,[],2) | samples < 0.999*min(avgs,[],2)));
            end
            % Unmark p > 1 cells if their extrema are smooth:
            for element = mesh.elements(isTroubled(:,:,priority) & [mesh.elements.dofCount] > 2)
                % Approximate 2nd and 3rd Legendre coefficients:
                coefs = element.getLegendre(2:3)./[1 3];
                if element.elementL.dofCount > 1
                    coefsL = element.elementL.getLegendre(2);
                else
                    coefsL = 0;
                end
                if element.elementR.dofCount > 1
                    coefsR = element.elementR.getLegendre(2);
                else
                    coefsR = 0;
                end
                element.isTroubled(:,:,priority) = any(abs(3*coefs(:,2) - Limiter_legendre.minmod(3*coefs(:,2),coefs(:,1)-coefsL,coefsR-coefs(:,1))) > 1e-10);
            end
        end
    end
end