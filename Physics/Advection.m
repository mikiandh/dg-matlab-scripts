classdef Advection < Physics
    properties
        advSpeed
        ghostStates
    end
    methods
        %% Constructor
        function advection = Advection(a,ghostStates)
            advection.equationCount = 1;
            advection.advSpeed = 1;
            if nargin > 0
                advection.advSpeed = a;
            end
            if nargin > 1
                advection.ghostStates = ghostStates;
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
            convection = this.advSpeed.*basis.gradientMatrix; % k_ij = v_j · c_ji
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
        %% Riemann solver (exact)
        function [flux,S] = riemannFlux(this,stateL,stateR)
            S = this.advSpeed;
            if S > 0
                flux = S*stateL;
            else
                flux = S*stateR;
            end
        end
        %% Boundary conditions
        function applyBoundaryConditions(this,mesh)
            if isempty(this.ghostStates)
                % Periodic boundary conditions:
                mesh.elements(end).riemannR = this.riemannFlux(...
                    mesh.elements(end).stateR,mesh.elements(1).stateL);
                mesh.elements(1).riemannL = - mesh.elements(end).riemannR;
                mesh.edges{end}.computeTimeDeltas(this.advSpeed);
                mesh.edges{1}.computeTimeDeltas(this.advSpeed);
            else
                % Transmissive boundary conditions:
                mesh.elements(end).riemannR = this.riemannFlux(mesh.elements(end).stateR,this.ghostStates(2));
                mesh.edges{end}.computeTimeDeltas(this.advSpeed);
                mesh.elements(1).riemannL = - this.riemannFlux(this.ghostStates(1),mesh.elements(1).stateL);
                mesh.edges{1}.computeTimeDeltas(this.advSpeed);
            end
        end
%         %% Periodic boundary condition
%         function applyBoundaryConditions(this,mesh)
%             mesh.elements(end).riemannR = this.riemannFlux(...
%                 mesh.elements(end).stateR,mesh.elements(1).stateL);
%             mesh.elements(1).riemannL = - mesh.elements(end).riemannR;
%             % Local time-step sizes:
%             mesh.edges{end}.computeTimeDeltas(this.advSpeed);
%             mesh.edges{1}.computeTimeDeltas(this.advSpeed);
%         end
        %% Outlet boundary condition
%         function applyBoundaryConditionsWeakly(this,mesh)
%             mesh.elements(end).riemannR = this.flux(mesh.elements(end).stateR);
%             mesh.elements(1).riemannL = 0;
%             % Local time-step sizes:
%             mesh.edges{end}.computeTimeDeltas(this.advSpeed);
%             mesh.edges{1}.computeTimeDeltas(this.advSpeed);
%         end
    end
end