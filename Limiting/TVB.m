classdef TVB < Limiter
    % Slope limiter by Cockburn and Shu (1998 and 2001). Proven to be TVB 
    % in the element-wise means (TVD, if M = 0). May reduce accuracy to 1st
    % order at smooth extrema (if M is not chosen properly).
    %
    properties (SetAccess = immutable)
        M % user-definable sensitivity parameter (see Cockburn & Shu, 1998)
    end
    methods
        %% Constructor
        function this = TVB(varargin)
            % Superclass constructor:
            this@Limiter(varargin{:});
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'M',0);
            % Parse the inputs:
            parse(p,varargin{:});
            this.M = p.Results.M;
        end
        %% Apply (extension)
        function applyStage(this,mesh,solver)
            % Default limiting:
            applyStage@Limiter(this,mesh,solver);
            % Retrieve troubled elements:
            elements = mesh.elements([mesh.elements.isTroubled]);
            % Apply on every remaining element:
            for element = elements
                this.applyOnElement(element);
            end
        end
        %% Information (extension)
        function info = getInfo(this)
            info = getInfo@Limiter(this);
            info = sprintf('%s, M = %g',info,this.M);
        end
    end
    methods (Access = protected)
        %% Limit an element
        function applyOnElement(this,element)
            % Applies TVB limiting (based on the TVD limiter of Osher) to
            % limit the slope (and, if p > 1, also the curvature) of the 
            % solution in this element (in local characteristic variables).
            %
            % Considers a cell to be troubled when EITHER of its edge
            % intercepts is beyond the range allowed by its neighbor's 
            % averages (i.e. NOT a la 1998, but a la 2001 and later).
            %
            % For p > 1, the limited solution is the parabola with safest
            % values at each edge, as determined in the trouble detection 
            % step.
            %
            % Aliases:
            tol = this.M*element.dx^2;
            basis = element.basis;
            % Get current cell Legendre coefficients:
            coefs = basis.getLegendre(element);
            % Cell slopes (via finite difference approximations):
            element.interpolateStateAtEdges;
            v1L = coefs(:,1) - element.stateL; % unsafe, left-sided
            v1R = element.stateR - coefs(:,1); % unsafe, right-sided
            u1L = coefs(:,1) - basis.getLegendre(element.elementL,1); % safe, left-sided
            u1R = basis.getLegendre(element.elementR,1) - coefs(:,1); % safe, right-sided
            % Compute the mapping to/from local characteristic variables:
            [~,L,R] = this.physics.getEigensystemAt(coefs(:,1));
            % Retrieve limited slopes:
            u1L = L*u1L; u1R = L*u1R;
            w1L = R*this.minmod(L*v1L,u1L,u1R,tol); % left-sided
            w1R = R*this.minmod(L*v1R,u1L,u1R,tol); % right-sided
            % Determine troubled rows and columns:
            rows = find(abs(w1L - v1L) > 1e-10 | abs(w1R - v1R) > 1e-10);
            element.isLimited(rows,2:end) = true; % flag all limited DOFs as such
            % Enforce safe edge intercepts in this element's solution:
            aux = zeros(size(coefs)); % preallocate to zero
            aux(rows,1) = .5*(w1R(rows) + w1L(rows)); % safe 2nd Legendre coefficient
            aux(rows,2) = 0*.5*(w1R(rows) - w1L(rows)); % safe 3rd Legendre coefficient (see A5-XI, 02/03/2020)
            coefs(:,2:end) = aux(:,1:basis.basisCount-1); % overwrite limited Legendre coefficients (only)
            basis.setLegendre(element,coefs(rows,:),rows); % overwrite all DOFs of limited state vector components
        end
    end
    methods (Static)
        %% Modified minmod operator (Shu, 1987)
        function A = minmod(A,B,C,tol)
            % Modified minmod function for 3 matrices and a given tolerance
            % (tol = M*dx^2). All arrays must be of the same sizes. Returns
            % a 2D array. "Inlined" for speed.
            %
            A = (abs(A) <= tol).*A + (abs(A) > tol).*(sign(A) == sign(B) & sign(A) == sign(C)).*sign(A).*min(abs(A),min(abs(B),abs(C)));
        end
    end
end