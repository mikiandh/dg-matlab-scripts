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
    end
    methods
        %% Sensor
        function apply(this,mesh,solver)
            % Apply default sensor first:
            apply@Sensor(this,mesh)
            % Evaluate state at all edges:
            mesh.elements.interpolateStateAtEdges
            % Check all possibly troubled elements:
            for element = mesh.elements([mesh.elements.isTroubled])
                % Compute jump in indicator variable across inflow edges:
                if solver.physics.getVelocityAt(element.stateL) > 0 % left edge inflow
                    I = abs(element.stateL(1) - element.elementL.stateR(1));
                else
                    I = 0; % reset
                end
                if solver.physics.getVelocityAt(element.stateR) < 0 % right edge inflow
                    I = I + abs(element.stateR(1) - element.elementR.stateL(1));
                end
                % Normalize it by an estimated baseline convergence rate:
                I = I/(abs(element.getLegendre(1,1))*element.dx^(.5*element.dofCount));
                % Set sensor status:
                element.isTroubled = I > this.treshold;
            end
        end
    end
end