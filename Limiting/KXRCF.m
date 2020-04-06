classdef KXRCF < Sensor
    %
    % Sensor by Krivodonova et al., 2004. Based on the superconvergence 
    % property that DG has at outflow boundaries of cells where the exact
    % solution is smooth, where order of accuracy rises to 2(p+1).
    %
    % Employs the first state vector component as indicator variable.
    % Velocity-like variable/quantity is physics-dependent.
    %
    properties (Constant)
        treshold = 1 % activation treshold
        indicator = struct('Euler',[1 3]) % indicator variables override
    end
    methods
        %% Sensor
        function apply(this,mesh,solver)
            % Apply default sensor first:
            apply@Sensor(this,mesh)
            % Evaluate state at all edges:
            mesh.elements.interpolateStateAtEdges
            % Set components of indicator variables:
            try
                i = this.indicator.(class(solver.physics));
            catch
                i = 1; % default (first component)
            end
            % Check all possibly troubled elements:
            for element = mesh.elements([mesh.elements.isTroubled])
                % Compute jump in indicator variable across inflow edges:
                if solver.physics.getVelocityAt(element.stateL) > 0 % left edge inflow
                    I = abs(element.stateL(i) - element.elementL.stateR(i));
                else
                    I = 0; % reset
                end
                if solver.physics.getVelocityAt(element.stateR) < 0 % right edge inflow
                    I = I + abs(element.stateR(i) - element.elementR.stateL(i));
                end
                % Normalize it by an estimated baseline convergence rate:
                I = I./(abs(element.getLegendre(1,i))*element.dx^(.5*element.dofCount));
                % Set sensor status:
                element.isTroubled = any(I > this.treshold);
            end
        end
    end
end