classdef Burgers < Physics
    properties (Constant, Hidden)
        equationCount = 1
    end
    properties
        viscosity = 0
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
        %% Jacobian eigen-decomposition
        function [D,L,R] = getEigensystemAt(state1,state2)
            % Returns the eigenvalue and eigenvector matrices evaluated at,
            % either:
            %
            % A) the given state vector
            % B) the "generalized Roe average" between two given states
            %
            if nargin == 1
                D = state1;
            else
                D = .5*(state1+state2); % Roe average for Burgers' reduces to arithmetic average (see Toro 2009, section 11.2)
            end
            L = 1;
            R = 1;
        end
        %% Advection velocity
        function v = getVelocityAt(state)
            v = .5*state;
        end
    end
    methods
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
    end
end