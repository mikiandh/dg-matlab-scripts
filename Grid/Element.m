classdef Element < handle
    % Class that represents the fundamental geometrical unit in which the
    % solution lives: the (DG) element, a.k.a patch.
    properties
        % DG:
        x
        xL
        xR
        dx
        basis % handle to a basis sub-class
        states % state vector at each degree of freedom (node or mode)
        fluxes % flux vector at each degree of freedom (Fletcher's group formulation)
        stateL % state vector at left edge
        stateR
        fluxL % flux vector at left edge
        fluxR
        riemannL % Riemann (upwind) flux vector at left edge
        riemannR
        edgeL % handle to its left edge
        edgeR
        elementR % handle to the closest element to its right
        elementL
        residuals % residual (i.e. time derivative) vector at nodes
        dofCount % number of degrees of freedom per element, per equation
        % Limiting:
        isTroubled = false % to be ignored by any limiter if set to false (set by the sensor)
        isLimited % true if modified by a limiter (row: system comp.; column: basis comp.)
        diffusions % artificial diffusion operator (cell array)
        antidiffusiveFluxes % sparse 2D array of FCT antidiffusive flux vectors along pairs of modes (with zero padding)
        maxima % 2D array of distances to local maximum state values, column: control point; row: system component
        minima % idem, for minima
        % SSP RK:
        extraStates
        extraResiduals
        localTimeDelta
    end
    methods
        %% Constructor (vector)
        function these = Element(bases,x)
            if nargin > 0
                validateattributes(x,"numeric",{'numel',numel(bases)+1})
                if numel(bases) > 1
                    these(1,numel(bases)) = Element;
                end
                [these.dofCount] = bases.basisCount;
                aux = num2cell(bases); [these.basis] = aux{:};
                aux = num2cell(x(1:end-1)); [these.xL] = aux{:};
                aux = num2cell(x(2:end)); [these.xR] = aux{:};
                aux = num2cell(mean([x(1:end-1);x(2:end)])); [these.x] = aux{:};
                aux = num2cell(diff(x)); [these.dx] = aux{:};
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
        %% Get DOF coordinates of this element (in physical space)
        function x = getDofCoords(this)
            x = this.basis.dofCoords;
            x = this.mapFromReference(x);
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
        %% Correct element states by adding antidiffusive fluxes
        function applyAntidiffusiveFluxes(this)
            % Updates the state vectors on an element by adding its
            % antidiffusive fluxes, weighted appropriately by the lumped
            % mass matrix entries.
            %
            % Aliases:
            N = this.basis.basisCount;
            masses = this.basis.lumpedMassMatrixDiagonal;
            % Loop over equations:
            for i = 1:size(this.states,1)
                f = reshape(this.antidiffusiveFluxes(i,:),N,N);
                this.states(i,:) = this.states(i,:) + sum(f,2)'./masses;
            end
        end
        %% Remove antidiffusive fluxes from element states
        function removeAntidiffusiveFluxes(this,betas)
            % Updates the state vectors on an element by removing a faction
            % of its antidiffusive fluxes, weighted appropriately by the 
            % lumped mass matrix entries.
            %
            % Aliases:
            N = this.basis.basisCount;
            masses = this.basis.lumpedMassMatrixDiagonal;
            % Loop over equations:
            for i = 1:size(this.states,1)
                f = betas.*reshape(this.antidiffusiveFluxes(i,:),N,N);
                this.states(i,:) = this.states(i,:) - sum(f,2)'./masses;
            end
        end
        %% Get Legendre coefficients of this element
        function Q = getLegendre(this,varargin)
            Q = this.basis.getLegendre(this,varargin{:});
        end
        %% Set Legendre coefficients of this element
        function setLegendre(this,varargin)
            this.basis.setLegendre(this,varargin{:});
        end
    end
end