classdef Burgers < Physics
    properties
        viscosity
        ghostStates
    end
    methods (Static)
        %% Convection operator
        function convection = getConvectionMatrix(states,basis)
            % Returns the discrete convection operator to be applied to the
            % modal coefficients of the discretizated solution.
            %
            % Arguments
            %  states: row array of basis function coefficients
            %  basis: finite-dimensional discretization basis
            % Output
            %  convection: discrete convection matrix, such that k_ij = v_j · c_ji
            convection = .5*states'.*basis.gradientMatrix;
            % convection = .5*spdiags(states',0,basis.basisCount,basis.basisCount)*basis.gradientMatrix; % maybe faster?
        end
        %% Flux function
        function flux = flux(state)
            flux = 0.5*state.^2;
        end
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
        %% Flux jacobian
        function A = jacobian(this,element)
            % Returns the Jacobian (in this case, a scalar) of the
            % conserved fluxes.
            A = this.advSpeed;
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
        %% Boundary conditions
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