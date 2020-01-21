classdef Burgers < Physics
    properties
        viscosity
        ghostStates
    end
    methods
        %% Constructor
        function burgers = Burgers(ghostStates,epsilon)
            if nargin < 2
                epsilon = 0;
            end
            if nargin > 0
                burgers.ghostStates = ghostStates;
            end
            burgers.equationCount = 1;
            burgers.viscosity = epsilon;
        end
        %% Jacobian matrix
        function A = getJacobian(~,q)
            % Returns the Jacobian matrices corresponding to given state
            % vectors.
            %
            % Argument
            %  q: 2D array of column state vectors
            % Output
            %  A: 3D array of Jacobian matrices (1x1), evaluated at each 
            %     state vector (page index)
            A = reshape(q,1,1,[]);
        end
        %% Flux function
        function flux = flux(~,state)
            flux = 0.5*state.^2;
        end
        %% Riemann solver (exact, with entropy fix)
        function [flux,S] = riemannFlux(this,stateL,stateR)
            S = [stateL stateR]';
            if stateL > 0 
                if stateR >= 0 % shock or expansion towards the right
                    flux = this.flux(stateL);
                elseif stateL > -stateR % shock towards the right
                    flux = this.flux(stateL);
                else % shock towards the left
                    flux = this.flux(stateR);
                end
            else
                if stateR < 0 % shock or expansion towards the left
                    flux = this.flux(stateR);
                else % transonic expansion (entropy fix)
                    flux = 0;
                end
            end
        end
        function applyBoundaryConditions(this,mesh)
            if isempty(this.ghostStates)
                % Periodic boundary conditions:
                [mesh.elements(end).riemannR,waveSpeeds] = this.riemannFlux(...
                    mesh.elements(end).stateR,mesh.elements(1).stateL);
                mesh.elements(1).riemannL = - mesh.elements(end).riemannR;
                mesh.edges{end}.computeTimeDeltas(waveSpeeds);
                mesh.edges{1}.computeTimeDeltas(waveSpeeds);
            else
                % Transmissive boundary conditions:
                [mesh.elements(end).riemannR,waveSpeeds] = this.riemannFlux(mesh.elements(end).stateR,this.ghostStates(2));
                mesh.edges{end}.computeTimeDeltas(waveSpeeds);
                [mesh.elements(1).riemannL,waveSpeeds] = this.riemannFlux(this.ghostStates(1),mesh.elements(1).stateL);
                mesh.elements(1).riemannL = - mesh.elements(1).riemannL;
                mesh.edges{1}.computeTimeDeltas(waveSpeeds);
            end
        end
    end
end