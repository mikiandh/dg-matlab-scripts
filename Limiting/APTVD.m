classdef APTVD < Sensor
    %
    % Accuracy-preserving TVD sensor based on Wang, 2009. Identifies local 
    % smoothness in the first derivatives of the solution (indicating a
    % smooth extrema). Implemented using Legendre coefficients in lieu of
    % cell-averaged derivatives (motivated by Krivodonova 2007).
    %
    methods
        %% Sensor
        function apply(this,mesh,~)
            % Apply default sensor first:
            apply@Sensor(this,mesh);
            % Mark all cells containing one or more local extrema:
            for element = mesh.elements([mesh.elements.isTroubled])
                % Get min/max cell average array:
                avgs = [element.basis.getLegendre(element,1)...
                    element.elementL.basis.getLegendre(element.elementL,1)...
                    element.elementR.basis.getLegendre(element.elementR,1)];
                % Get solution samples:
                samples = element.interpolateStateAtCoords([element.xL element.getDofCoords' element.xR]);
                % Unmark if no local extremum is present:
                element.isTroubled = any(any(samples > 1.001*max(avgs,[],2) | samples < 0.999*min(avgs,[],2)));
            end
            % Unmark p > 1 cells if their extrema are smooth:
            for element = mesh.elements([mesh.elements.isTroubled] & [mesh.elements.dofCount] > 2)
                % Approximate 2nd and 3rd Legendre coefficients:
                coefs = element.basis.getLegendre(element,2:3)./[1 3];
                if element.elementL.dofCount > 1
                    coefsL = element.elementL.basis.getLegendre(element.elementL,2);
                else
                    coefsL = 0;
                end
                if element.elementR.dofCount > 1
                    coefsR = element.elementR.basis.getLegendre(element.elementR,2);
                else
                    coefsR = 0;
                end
                element.isTroubled = any(abs(3*coefs(:,2) - BDF.minmod(3*coefs(:,2),coefs(:,1)-coefsL,coefsR-coefs(:,1))) > 1e-10);
            end
        end
    end
end