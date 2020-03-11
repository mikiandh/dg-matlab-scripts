classdef APTVD < Sensor
    %
    % Accuracy-preserving TVD sensor based on Wang, 2009. Identifies local 
    % smoothness in the first derivatives of the solution (indicating a
    % smooth extrema). Implemented using Legendre coefficients in lieu of
    % cell-averaged derivatives (motivated by Krivodonova 2007).
    %
    properties
        beta
    end
    methods
        %% Constructor
        function this = APTVD(varargin)
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'beta',1.5,@(x) x >= 1 && x <= 2);
            % Parse the inputs:
            parse(p,varargin{:});
            this.beta = p.Results.beta;
        end
        %% Sensor
        function apply(this,mesh,~)
            % Apply default sensor first:
            apply@Sensor(this,mesh);
            % Mark all cells containing one or more local extrema:
            for element = findobj(mesh.elements,'isTroubled',true)'
                % Get min/max cell average array:
                avgs = [element.basis.getLegendre(element,1)...
                    element.edgeL.elementL.basis.getLegendre(element.edgeL.elementL,1)...
                    element.edgeR.elementR.basis.getLegendre(element.edgeR.elementR,1)];
                % Get solution samples at edges:
                element.interpolateStateAtEdges;
                samples = [element.stateL element.stateR];
                % Unmark if no local extremum is present:
                element.isTroubled = any(any(samples > 1.001*max(avgs,[],2) | samples < 0.999*min(avgs,[],2)));
            end
            % Unmark p > 1 cells if their extrema are smooth:
            for element = findobj(mesh.elements,'isTroubled',true,'-not','dofCount',2)'
                % Aliases:
                elementL = element.edgeL.elementL;
                elementR = element.edgeR.elementR;
                % Approximate 1st and 2nd cell-averaged derivatives:
                ders = element.basis.getLegendre(element,2:3)./[1 3];
                if elementL.dofCount > 1
                    dersL = elementL.basis.getLegendre(elementL,2);
                else
                    dersL = 0;
                end
                if elementR.dofCount > 1
                    dersR = elementR.basis.getLegendre(elementR,2);
                else
                    dersR = 0;
                end
                % Compare 2nd derivative with 1st derivative differences:
                element.isTroubled = any(ders(:,2) ~= BDF.minmod(ders(:,2),...
                    this.beta*element.dx^2/(element.dx + elementL.dx)*(ders(:,1)/element.dx - dersL/elementL.dx),...
                    this.beta*element.dx^2/(element.dx + elementR.dx)*(dersR/elementR.dx - ders(:,1)/element.dx)));
            end
        end
    end
end