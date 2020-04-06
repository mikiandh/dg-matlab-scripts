classdef Advection < Physics
    properties (Constant, Hidden)
        equationCount = 1
        controlVars = 1
    end
    properties (SetAccess = immutable)
        advSpeed = 1
    end
    methods
        %% Constructor
        function this = Advection(a)
            if nargin
                validateattributes(a,"numeric",{'finite'})
                this.advSpeed = a;
            end
        end
        %% Convection operator
        function convection = getConvectionMatrix(this,~,basis)
            % Returns the discrete convection operator to be applied to the
            % modal coefficients of the discretizated solution.
            %
            % Arguments
            %  states: (unused) row array of basis function coefficients
            %  basis: finite-dimensional discretization basis
            % Output
            %  convection: discrete convection operator
            convection = this.advSpeed.*basis.gradientMatrix; % k_ij = v_j � c_ji
        end
        %% Flux function
        function flux = flux(this,state)
            flux = this.advSpeed*state;
        end
        %% Flux jacobian
        function A = getJacobian(this,~,~)
            % Returns the "edge Jacobian" (in this case, a scalar) of the
            % conserved fluxes, evaluated at the "virtual edge" connecting 
            % the 1st given state with 2nd.
            %
            % Arguments
            %  statesL: (not used) first state
            %  statesR: (not used) second state
            %
            A = this.advSpeed;
        end
        %% Jacobian eigen-decomposition
        function [D,L,R] = getEigensystemAt(this,~,~)
            % Returns the eigenvalue and eigenvector matrices evaluated at,
            % either:
            %
            % A) the given state vector
            % B) the "generalized Roe average" between two given states
            %
            D = this.advSpeed;
            L = 1;
            R = 1;
        end
        %% Riemann solver (exact)
        function [flux,S] = riemannFlux(this,stateL,stateR)
            S = this.advSpeed;
            if S > 0
                flux = S*stateL;
            else
                flux = S*stateR;
            end
        end
        %% Advection velocity
        function v = getVelocityAt(this,~)
            v = this.advSpeed;
        end
        %% Information (extension)
        function info = getInfo(this)
            info = sprintf('%s, a = %g',this.getInfo@Physics,this.advSpeed);
        end
    end
end