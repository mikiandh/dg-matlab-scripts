classdef Element < handle
    properties
        xL
        xR
        dx
        basis % handle to a basis sub-class
        states % state vector at each degree of freedom (node or mode)
        fluxes % flux vector at nodes (if any)
        stateL % state vector at left edge
        stateR
        fluxL % flux vector at left edge
        fluxR
        riemannL % Riemann (upwind) flux vector at left edge
        riemannR
        residuals % residual (i.e. time derivative) vector at nodes
        extraStates
        extraResiduals
        localTimeDelta
        limiterHistory
        dofCount % number of degrees of freedom per element, per equation
    end
    methods
        %% Constructor
        function element = Element(x1,x2,basis)
            if nargin > 0
                element.xL = x1;
                element.xR = x2;
                element.dx = x2-x1;
                element.basis = basis;
                element.dofCount = basis.basisCount;
            end
        end
        %% Affine map (from physical to reference element space)
        function [xi,jac] = mapToReference(this,x)
            jac = 2/this.dx;
            xi = jac*(x-this.xL) - 1;
        end
        %% Affine map (from reference element space to physical)
        function [x,jac] = mapFromReference(this,xi)
            jac = 0.5*this.dx;
            x = jac*(xi+1) + this.xL;
        end
        %% Get node coordinates of this element (in physical space)
        function x = getNodeCoords(this)
            x = this.basis.nodeCoords;
            x = this.mapFromReference(x);
        end
        %% Get quadrature point coordinates of this element (in physical space)
        function x = getGaussCoords(this)
            x = this.basis.getGaussCoords;
            x = this.mapFromReference(x);
        end
        %% Get breakpoint coordinates within this element/patch (in physical space)
        function x = getBreakCoords(this)
            x = this.basis.breakCoords;
            x = this.mapFromReference(x);
        end
        %% Get control coordinates within this element/patch (in physical space)
        function x = getControlCoords(this)
            x = this.basis.controlCoords;
            x = this.mapFromReference(x);
        end
        %% Convert element from Lagrange (nodal) to Legendre (modal) basis
        function nodalToModal(this)
            this.states = this.states*this.basis.invVandermonde;
        end
        %% Convert element from Legendre (modal) to Lagrange (nodal) basis
        function modalToNodal(this)
            this.states = this.states*this.basis.vandermonde;
        end
        %% Interpolate state at coordinates
        function q = interpolateStateAtCoords(this,x)
            % Argument
            %  x: M locations in physical space where to sample this 
            %     element's fin. dim. approx. to the solution (row array).
            % Return
            %  q: state vector evaluated at the M sample locations 'x' 
            %     (2D matrix, nEqs x M).
            %
            x = this.mapToReference(x);
            q = this.states*this.basis.sampleAt(x);
        end
        %% Interpolate flux at coordinates
        function f = interpolateFluxAtCoords(this,x)
            % Argument
            %  x: M locations in physical space where to sample this 
            %     element's fin. dim. approx. to the solution (row array).
            % Return
            %  f: flux vector evaluated at the M sample locations 'x' 
            %     (2D matrix, nEqs x M).
            %
            x = this.mapToReference(x);
            f = this.fluxes*this.basis.sampleAt(x);
        end
        %% Interpolate residual at coordinates
        function q = interpolateResidualAtCoords(this,x)
            % Argument
            %  x: M locations in physical space where to sample this 
            %     element's fin. dim. approx. to the solution (row array).
            % Return
            %  q: residual vector evaluated at the M sample locations 'x' 
            %     (2D matrix, nEqs x M).
            %
            % Will produce an error (matrix dimension mismatch) if 
            % residudals have not been initialized prior to this call.
            %
            x = this.mapToReference(x);
            q = this.residuals*this.basis.sampleAt(x);
        end
        %% Interpolate state at element edges
        function interpolateStateAtEdges(this)
            % Evaluates and stores the state vector at the element's edges.
            this.stateL = this.states*this.basis.left;
            this.stateR = this.states*this.basis.right;
        end
        %% Interpolate flux at element edges
        function interpolateFluxAtEdges(this)
            % Evaluates and stores the flux vector at the element's edges.
            this.fluxL = this.fluxes*this.basis.left;
            this.fluxR = this.fluxes*this.basis.right;
        end
        %% Apply flux function to states
        function computeFluxesFromStates(this,physics)
            this.fluxes = physics.flux(this.states);
        end
        %% Update element residuals according to its basis
        function computeResiduals(this,physics)
            this.basis.computeResiduals(this,physics);
        end
    end
end