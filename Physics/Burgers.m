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
        %% Solution breakdown
        function [t,x,x0] = getBreakingPoint(f,L,varargin)
            % Given an initial condition function, tries to find the
            % earliest breaking point in space-time, for t > 0.
            %
            % Inputs
            %  f: an initial condition function (of one variable)
            %  x: left and right boundary positions
            %  isPeriodic: assume that f is periodic outside [x1,x2] or not
            %  x0: custom array of seed points
            %
            % Outputs
            %  t: breaking time
            %  x: one breaking location (of possibly many)
            %  x0: origin of 2 (of possibly more) charac. meeting at (t,x)
            %
            % Parser:
            validateattributes(f,{'function_handle'},{'nonempty'});
            p = inputParser;
            %%%addParameter(p,'x0',linspace(L(1),L(end),10))
            addParameter(p,'x0',L(1) + diff(L)*[0 rand(1,8) 1])
            addParameter(p,'isPeriodic',false)
            parse(p,varargin{:});
            % Periodic extension, if requested:
            if p.Results.isPeriodic
                f = @(x) f(L(1) + mod(x - L(1),L(end) - L(1)));
            end
            % Line-line intersections:
            function t = breakdown(x0)
                % Given an array of starting points, returns the earliest
                % intersection time, i.e. the breaking time of the
                % solution. Discards negative intersection times, parallel
                % lines, and self-intersections. Assumes periodic initial
                % condition.
                %
                % I/O
                %  x0: 1D array of origin points
                %  t: breakdown time
                t = (x0.' - x0)./(f(x0) - f(x0).');
                t(~(t > 0)) = [];
                t = min([t(t > 0) inf]); % always returns something
            end
            [x0,t] = fminsearch(@breakdown,p.Results.x0); % breaking time
            x = f(x0).*t + x0;
            x0 = x0(find(sum(x == x') > 1,2));
            x = x(find(sum(x == x') > 1,1));
        end
        %% Method of characteristics (robust but slow)
        function X = MOC_robust(t,x,f,L)
            % Solves the given function via the M.O.C.; does NOT try to
            % detect breaking points: coming across one will result in
            % undefined behaviour.
            %
            % Slow but robust version (crashes more gracefully near
            % breaking points than 'MOC').
            %
            % Periodic extension, if interval bounds were given:
            if nargin > 3
                f = @(x) f(L(1) + mod(x - L(1),L(end) - L(1)));
            end
            % Try to recover the origin point of each cha. reaching (x,t):
            [X,T] = meshgrid(x,t);
            options = optimset('Display','off');
            for i = 1:numel(X) 
                X(i) = fzero(@(X0) X(i) - f(X0)*T(i) - X0,X(i),options);
            end
            X = f(X);
        end
        %% Method of characteristics (fast but delicate)
        function u = MOC(t,x,f,L)
            % Solves the given function via the M.O.C.; does NOT try to
            % detect breaking points: coming across one will result in
            % undefined behaviour.
            %
            % Only accepts a scalar time instant; usually, much
            % faster than 'MOC_robust'.
            %
            % Periodic extension, if interval bounds were given:
            if nargin > 3
                f = @(x) f(L(1) + mod(x - L(1),L(end) - L(1)));
            end
            % Fixed point iteration:
            u0 = f(x);
            for iter = 1:1000
                u = f(x - t.*u0);
                if max(abs(u - u0)) < 1e-6
                    return
                else
                    u0 = u;
                end
            end
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