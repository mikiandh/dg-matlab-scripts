classdef TVB < Limiter
    %
    % Slope limiter as reported in Cockburn and Shu, 2001. Proven to be TVB
    % in the element-wise means (TVD if M = 0). Known to cause loss of 
    % accuracy at smooth extrema unless M is set appropriately.
    % 
    % For p > 1, the solution is L2-projected on a linear basis, and 
    % the resulting slope is used as in the linear definition of this
    % limiter (i.e. the 2nd expansion coefficient of the Legendre basis is
    % always used as the unlimited slope, regardless of p).
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
            elements = findobj(mesh.elements,'isTroubled',true);
            % Discard any element with a non-modal basis:
            elements(arrayfun(@(x) ~x.basis.isModal,elements)) = [];
            % Apply on every remaining element:
            for element = elements'
                this.applyOnElement(element);
            end
        end
    end
    methods (Access = protected)
        function applyOnElement(this,element)
            % Interpolate solution at element edges:
            element.interpolateStateAtEdges;
            % Aliases:
            uL = element.edgeL.elementL.states(:,1); % left cell means
            uR = element.edgeR.elementR.states(:,1); % right cell means
            u = element.states(:,1); % current cell means
            vL = element.stateL; % unlimited left edge intercept
            vR = element.stateR; % unlimited right edge intercept
            % Limited edge intercepts:
            wL = u - this.minmod(u-vL,u-uL,uR-u); % left edge
            wR = u + this.minmod(vR-u,u-uL,uR-u); % right edge
            % Determine troubled components:
            eqs = find(wL ~= vL | wR ~= vR);
            % Overwrite troubled slopes with limited slopes:
            for i = eqs'
                s = .5*(vR(i) - vL(i)); % unlimited slope
                sL = u(i) - uL(i); % left-limited slope
                sR = uR(i) - u(i); % right-limited slope
                element.states(i,2) = this.minmod(s,sL,sR); % most conservative slope
                element.states(i,3:end) = 0; % reduce to linear basis
                element.isLimited(i,2:element.dofCount) = true;
            end
        end
    end
    methods (Static, Access = protected)
        %% Minmod (3 column array arguments)
        function d = minmod(a,b,c)
            d = zeros(length(a),1);
            M = abs(horzcat(a,b,c));
            a = sign(a);
            ids = a == sign(b) & a == sign(c);
            d(ids) = min(M(ids,:),[],2);
            d = a.*d;
        end
    end
end