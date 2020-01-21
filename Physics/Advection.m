classdef Advection < Physics
    properties
        advSpeed
    end
    methods
        %% Constructor
        function advection = Advection(a)
            advection.equationCount = 1;
            advection.advSpeed = 1;
            if nargin > 0
                advection.advSpeed = a;
            end
        end
        %% Jacobian matrix
        function A = getJacobian(this,q)
            % Returns the Jacobian matrices corresponding to given state
            % vectors.
            %
            % Argument
            %  q: 2D array of column state vectors
            % Output
            %  A: 3D array of Jacobian matrices, evaluated at each state
            %  vector (page index)
            A = repmat(this.advSpeed,1,1,length(q));
        end
        %% Flux function
        function flux = flux(this,state)
            flux = this.advSpeed*state;
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
        %% Periodic boundary condition
        function applyBoundaryConditions(this,mesh)
            mesh.elements(end).riemannR = this.riemannFlux(...
                mesh.elements(end).stateR,mesh.elements(1).stateL);
            mesh.elements(1).riemannL = - mesh.elements(end).riemannR;
            % Local time-step sizes:
            mesh.edges{end}.computeTimeDeltas(this.advSpeed);
            mesh.edges{1}.computeTimeDeltas(this.advSpeed);
        end
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