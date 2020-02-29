classdef TVB < Limiter
    %
    % Slope limiter as reported in Cockburn and Shu, 1998. Proven to be TVB
    % in the element-wise means (TVD, if M = 0). May reduce accuracy to 1st
    % order at smooth extrema (if M is not chosen properly).
    %
    properties
        M % user-definable sensitivity parameter (see Cockburn & Shu, 2001)
    end
    properties (Access = protected)
        physics % necessary for the characteristic projection
    end
    methods
        %% Constructor
        function this = TVB(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'M',0);
            % Parse the M parameter:
            parse(p,varargin{:});
            this.M = p.Results.M;
        end
        %% Apply
        function apply(this,mesh,~)
            % Apply default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
            % Project troubled elements and neighbors
            % Apply on every remaining element:
            for element = elements
                this.applyOnElement(element);
            end
        end
        %% Information
        function info = getInfo(this)
            info = sprintf('TVB, M = %g',this.M);
        end
    end
    methods (Access = protected)
        function applyOnElement(this,element)
            % Applies TVB limiting (based on the TVD limiter of Osher) to
            % limit the gradient of the solution in this element (in local
            % characteristic variables).
            %
            % Aliases:
            tol = this.M*element.dx^2;
            basis = element.basis;
            % Cell averages:
            v0 = basis.getLegendre(element,1);
            % Cell slopes (finite difference approximations):
            element.interpolateStateAtEdges;
            v1L = v0 - element.stateL; % unsafe, left-sided
            v1R = element.stateR - v0; % unsafe, right-sided
            u1L = v0 - basis.getLegendre(element.edgeL.elementL,1); % safe, left-sided
            u1R = basis.getLegendre(element.edgeR.elementR,1) - v0; % safe, right-sided
            % Retrieve limited slopes:
            w1L = this.modminmod([v1L u1L u1R],tol); % left-sided
            w1R = this.modminmod([v1R u1L u1R],tol); % right-sided
            % Determine troubled rows and columns:
            rows = find(w1L ~= v1L | w1R ~= v1R);
            cols = 2:basis.basisCount;
            % Overwrite troubled entries with limited ones:
            vals = zeros(length(rows),basis.basisCount+1);
            vals(:,2) = .5*(w1R(rows) + w1L(rows)); % safe slopes 
            vals(:,3) = .5*(w1R(rows) - w1L(rows)); % safe curvatures
            basis.setLegendre(element,vals(:,cols),cols,rows);
            element.isLimited(rows,cols) = true; % flag all limited modes
        end
    end
    methods (Static)
        %% Modified minmod operator (Shu, 1987)
        function outs = modminmod(args,tol)
            % Modified minmod function for a 2D array input of 3 columns
            % and a given tolerance (tol = M*dx^2). Returns a column array.
            %
            % Pre-process inputs:
            signs = sign(args);
            args = abs(args);
            % Curvature check:
            PASS = args(:,1) <= tol;
            % Slope check:
            FAIL = signs(:,1) ~= signs(:,2) | signs(:,1) ~= signs(:,3);
            % Minmod operator: 
            outs = signs(:,1).*min(args,[],2);
            outs(FAIL) = 0;
            outs(PASS) = signs(PASS,1).*args(PASS,1);
        end
    end
end