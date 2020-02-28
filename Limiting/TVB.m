classdef TVB < Limiter
    %
    % Slope limiter as reported in Cockburn and Shu, 2001. Proven to be TVB
    % in the element-wise means (TVD, if M = 0). May reduce accuracy to 1st
    % order at smooth extrema (if M is not chosen properly).
    %
    properties
        M % user-definable sensitivity parameter (see Cockburn & Shu, 2001)
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
            % Applies TVB limiting (based on the TVD limiter of Osher) a la
            % Cockburn & Shu, 2001.
            %
            % Minmod tolerance:
            tol = this.M*element.dx^2;
            % Cell averages:
            uL = element.basis.getLegendre(element.edgeL.elementL,1); % left cell
            uR = element.basis.getLegendre(element.edgeR.elementR,1); % right cell
            u = element.basis.getLegendre(element,1); % current cell
            % Unlimited edge intercepts:
            element.interpolateStateAtEdges;
            vL = element.stateL; % left edge
            vR = element.stateR; % right edge
            % Limited edge intercepts:
            wL = u - this.modminmod([u-vL u-uL uR-u],tol); % left edge
            wR = u + this.modminmod([vR-u u-uL uR-u],tol); % right edge
            % Determine troubled components:
            rows = find(wL ~= vL | wR ~= vR);
            cols = 3:element.dofCount; % higher order modes (above linear)
            % Overwrite troubled slopes with limited ones:
            for i = rows'
                s = .5*(vR(i) - vL(i)); % unlimited slope
                sL = u(i) - uL(i); % left-limited slope
                sR = uR(i) - u(i); % right-limited slope
                s = this.modminmod([s sL sR],tol); % most conservative slope
                element.basis.setLegendre(element,s,2,i); % override slope
                element.basis.setLegendre(element,0,cols,i); % truncate to linear
                element.isLimited(i,2:end) = true; % flag all limited modes
            end
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