classdef KXRCF < Sensor
    %
    % Sensro by Krivodonova et al., 2004. Based on a superconvergence 
    % property of DG at element edges where the solution is smooth.
    %
    methods
        %% Sensor
        function apply(this,mesh,~)
            % Apply default sensor first:
            apply@Sensor(this,mesh);
            % 
            %%% TO DO %%%
        end
    end
end